"""Build solids.txt: refined unit-edge coordinates for all convex regular-faced solids.
Usage: python3 prep.py DATA_DIR OUT"""
import json, math, os, sys
import numpy as np

CHIRAL = {"snub-cube", "snub-dodecahedron", "gyroelongated-triangular-bicupola", "gyroelongated-square-bicupola",
          "gyroelongated-pentagonal-bicupola", "gyroelongated-pentagonal-cupolarotunda", "gyroelongated-pentagonal-birotunda"}

def prism(m):
    R = 0.5 / math.sin(math.pi / m)
    v = [(R * math.cos(2 * math.pi * i / m), R * math.sin(2 * math.pi * i / m), z) for z in (0.0, 1.0) for i in range(m)]
    faces = [list(range(m)), list(range(m, 2 * m))]
    for i in range(m):
        j = (i + 1) % m
        faces.append([i, j, m + j, m + i])
    return np.array(v), faces

def antiprism(m):
    R = 0.5 / math.sin(math.pi / m)
    h = math.sqrt(1 - 0.25 / math.cos(math.pi / (2 * m)) ** 2)
    v = [(R * math.cos(2 * math.pi * i / m), R * math.sin(2 * math.pi * i / m), 0.0) for i in range(m)]
    v += [(R * math.cos(2 * math.pi * (i + 0.5) / m), R * math.sin(2 * math.pi * (i + 0.5) / m), h) for i in range(m)]
    faces = [list(range(m)), list(range(m, 2 * m))]
    for i in range(m):
        j = (i + 1) % m
        faces.append([i, j, m + i])
        faces.append([m + i, m + j, j])
    return np.array(v), faces

def constraints(faces):
    cons = {}
    for f in faces:
        n = len(f)
        for a in range(n):
            for b in range(a + 1, n):
                i, j = sorted((f[a], f[b]))
                cons[(i, j)] = math.sin(math.pi * (b - a) / n) / math.sin(math.pi / n)
    return cons

def refine(v, faces):
    v = v.copy()
    cons = constraints(faces)
    keys = list(cons)
    for it in range(50):
        J = np.zeros((len(keys), v.size))
        r = np.zeros(len(keys))
        for row, (i, j) in enumerate(keys):
            d = v[i] - v[j]
            l = np.linalg.norm(d)
            r[row] = l - cons[(i, j)]
            J[row, 3 * i:3 * i + 3] = d / l
            J[row, 3 * j:3 * j + 3] = -d / l
        if np.abs(r).max() < 1e-14:
            break
        dx = np.linalg.lstsq(J, -r, rcond=1e-10)[0]
        v += dx.reshape(v.shape)
    return v, np.abs(r).max()

def orient(v, faces):
    v = v - v.mean(axis=0)
    out = []
    for f in faces:
        f = list(f)
        p = v[f]
        fc = p.mean(axis=0)
        n = np.zeros(3)
        for i in range(len(f)):
            n += np.cross(p[i] - fc, p[(i + 1) % len(f)] - fc)
        if np.dot(n, fc) < 0:
            f = f[::-1]
        out.append(f)
    return v, out

def main():
    data, outp = sys.argv[1], sys.argv[2]
    solids = []
    for fn in sorted(os.listdir(data)):
        d = json.load(open(os.path.join(data, fn)))
        name = fn[:-5]
        v = np.array(d["vertices"], dtype=float)
        e = d["edges"][0]
        v = v / np.linalg.norm(v[e[0]] - v[e[1]])
        solids.append((name, v, d["faces"]))
    have = {s[0] for s in solids}
    pn = {3: "triangular", 5: "pentagonal", 6: "hexagonal", 8: "octagonal", 10: "decagonal"}
    an = {4: "square", 5: "pentagonal", 6: "hexagonal", 8: "octagonal", 10: "decagonal"}
    for m, nm in pn.items():
        if nm + "-prism" not in have:
            solids.append((nm + "-prism",) + prism(m))
    for m, nm in an.items():
        if nm + "-antiprism" not in have:
            solids.append((nm + "-antiprism",) + antiprism(m))
    with open(outp, "w") as out:
        worst = 0
        for name, v, faces in solids:
            v, err = refine(v, faces)
            worst = max(worst, err)
            v, faces = orient(v, faces)
            variants = [(name, v)]
            if name in CHIRAL:
                m = v.copy(); m[:, 0] = -m[:, 0]
                variants.append((name + "~", m))
            for nm, vv in variants:
                vv, ff = orient(vv, faces)
                out.write("solid %s %d %d\n" % (nm, len(vv), len(ff)))
                for p in vv:
                    out.write("%.17g %.17g %.17g\n" % tuple(p))
                for f in ff:
                    out.write("%d %s\n" % (len(f), " ".join(map(str, f))))
        print("solids", len(solids), "worst residual", worst)

main()
