"""Render an OFF file (with its .pieces companion) to a small PNG. numpy only.

Usage: render.py file.off out.png [size]
A z-buffer rasteriser: faces are fan-triangulated, flat shaded, coloured by the solid they
belong to. Meant for notification thumbnails, not for publication."""
import struct, sys, zlib
import numpy as np

PALETTE = [(0.36, 0.58, 0.84), (0.85, 0.44, 0.58), (0.47, 0.72, 0.52), (0.95, 0.74, 0.32),
           (0.60, 0.50, 0.80), (0.40, 0.75, 0.78), (0.90, 0.55, 0.35), (0.65, 0.65, 0.65)]
BG = (0.04, 0.05, 0.07)


def load_off(path):
    t = open(path).read().split()
    nv, nf = int(t[1]), int(t[2])
    v = np.array(t[4:4 + 3 * nv], dtype=float).reshape(nv, 3)
    faces, i = [], 4 + 3 * nv
    for _ in range(nf):
        n = int(t[i]); faces.append([int(x) for x in t[i + 1:i + 1 + n]]); i += n + 1
    return v, faces


def write_png(path, img):
    h, w, _ = img.shape
    raw = b"".join(b"\x00" + img[y].tobytes() for y in range(h))
    def chunk(tag, data):
        return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", zlib.crc32(tag + data) & 0xffffffff)
    with open(path, "wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)) + chunk(b"IDAT", zlib.compress(raw, 6)) + chunk(b"IEND", b""))


def render(off, out, size=512):
    v, faces = load_off(off)
    try:
        names = open(off + ".pieces").read().split("\n")
    except OSError:
        names = [""] * len(faces)
    kinds = sorted(set(names))
    v = v - v.mean(axis=0)
    # a fixed three-quarter view
    a, b = np.radians(25), np.radians(-20)
    ry = np.array([[np.cos(a), 0, np.sin(a)], [0, 1, 0], [-np.sin(a), 0, np.cos(a)]])
    rx = np.array([[1, 0, 0], [0, np.cos(b), -np.sin(b)], [0, np.sin(b), np.cos(b)]])
    p = v @ (rx @ ry).T
    s = 0.46 * size / np.linalg.norm(p, axis=1).max()
    px, py, pz = size / 2 + s * p[:, 0], size / 2 - s * p[:, 1], p[:, 2]
    img = np.empty((size, size, 3)); img[:] = BG
    zb = np.full((size, size), -np.inf)
    light = np.array([0.35, 0.5, 0.8]); light /= np.linalg.norm(light)
    for f, name in zip(faces, names):
        n = np.cross(p[f[1]] - p[f[0]], p[f[2]] - p[f[1]])
        ln = np.linalg.norm(n)
        if ln < 1e-12: continue
        shade = 0.35 + 0.65 * abs(float(n @ light) / ln)
        col = np.array(PALETTE[kinds.index(name) % len(PALETTE)]) * shade
        for k in range(1, len(f) - 1):
            i0, i1, i2 = f[0], f[k], f[k + 1]
            x0, x1 = int(max(0, np.floor(min(px[i0], px[i1], px[i2])))), int(min(size - 1, np.ceil(max(px[i0], px[i1], px[i2]))))
            y0, y1 = int(max(0, np.floor(min(py[i0], py[i1], py[i2])))), int(min(size - 1, np.ceil(max(py[i0], py[i1], py[i2]))))
            if x1 < x0 or y1 < y0: continue
            X, Y = np.meshgrid(np.arange(x0, x1 + 1) + 0.5, np.arange(y0, y1 + 1) + 0.5)
            d = (py[i1] - py[i2]) * (px[i0] - px[i2]) + (px[i2] - px[i1]) * (py[i0] - py[i2])
            if abs(d) < 1e-12: continue
            w0 = ((py[i1] - py[i2]) * (X - px[i2]) + (px[i2] - px[i1]) * (Y - py[i2])) / d
            w1 = ((py[i2] - py[i0]) * (X - px[i2]) + (px[i0] - px[i2]) * (Y - py[i2])) / d
            w2 = 1 - w0 - w1
            z = w0 * pz[i0] + w1 * pz[i1] + w2 * pz[i2]
            sub = zb[y0:y1 + 1, x0:x1 + 1]
            m = (w0 >= -1e-9) & (w1 >= -1e-9) & (w2 >= -1e-9) & (z > sub)
            sub[m] = z[m]
            img[y0:y1 + 1, x0:x1 + 1][m] = col
    write_png(out, (np.clip(img, 0, 1) * 255).astype(np.uint8))


if __name__ == "__main__":
    render(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 512)
