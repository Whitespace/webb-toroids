"""Independent check of an OFF file: closed oriented manifold, regular unit faces, genus, (A), (Q).
Usage: python3 check.py file.off [...]"""
import sys, math
import numpy as np
from scipy.spatial import ConvexHull

def load(path):
    t = open(path).read().split()
    nv, nf = int(t[1]), int(t[2])
    v = np.array(t[4:4 + 3 * nv], dtype=float).reshape(nv, 3)
    faces, i = [], 4 + 3 * nv
    for _ in range(nf):
        n = int(t[i]); faces.append([int(x) for x in t[i + 1:i + 1 + n]]); i += n + 1
    return v, faces

def check(path):
    v, faces = load(path)
    de = {}
    bad = 0
    normals = []
    for fi, f in enumerate(faces):
        n = len(f)
        p = v[f]
        c = p.mean(axis=0)
        for i in range(n):
            a, b = f[i], f[(i + 1) % n]
            if (a, b) in de: bad += 1
            de[(a, b)] = fi
            if abs(np.linalg.norm(v[a] - v[b]) - 1) > 1e-6: bad += 1
        r = np.linalg.norm(p - c, axis=1)
        if r.max() - r.min() > 1e-6: bad += 1
        nn = np.cross(p[1] - p[0], p[2] - p[1]); nn /= np.linalg.norm(nn)
        if np.abs((p - c) @ nn).max() > 1e-6: bad += 1
        normals.append(nn)
    for (a, b) in de:
        if (b, a) not in de: bad += 1
    E = len(de) // 2
    chi = len(v) - E + len(faces)
    aviol = sum(1 for (a, b), f in de.items() if a < b and np.dot(normals[f], normals[de[(b, a)]]) > 1 - 1e-9)
    hull = ConvexHull(v)
    # merge coplanar simplices into faces, collect true hull edges
    planes = {}
    for s, eq in zip(hull.simplices, hull.equations):
        key = tuple(np.round(eq, 5))
        planes.setdefault(key, []).append(s)
    edges_missing = 0; nedges = 0; sizes = []; irregular = 0
    for key, tris in planes.items():
        cnt = {}
        for s in tris:
            for i in range(3):
                e = tuple(sorted((s[i], s[(i + 1) % 3])))
                cnt[e] = cnt.get(e, 0) + 1
        border = [e for e, c in cnt.items() if c == 1]
        # chain border into polygon and drop collinear vertices
        nxt = {}
        for a, b in border:
            nxt.setdefault(a, []).append(b); nxt.setdefault(b, []).append(a)
        start = border[0][0]; poly = [start]; prev = None; cur = start
        while True:
            cand = [x for x in nxt[cur] if x != prev]
            n2 = cand[0]
            if n2 == start: break
            poly.append(n2); prev, cur = cur, n2
        corners = []
        m = len(poly)
        for i in range(m):
            a, b, c = v[poly[i - 1]], v[poly[i]], v[poly[(i + 1) % m]]
            if np.linalg.norm(np.cross(b - a, c - b)) > 1e-7: corners.append(poly[i])
        sizes.append(len(corners))
        pc = v[corners]; r = np.linalg.norm(pc - pc.mean(axis=0), axis=1)
        if r.max() - r.min() > 1e-6: irregular += 1
        for i in range(len(corners)):
            a, b = corners[i], corners[(i + 1) % len(corners)]
            nedges += 1
            if (a, b) not in de and (b, a) not in de: edges_missing += 1
    # embeddedness: no vertex of the mesh strictly inside another face's triangle fan is checked elsewhere (Rust SAT)
    hist = {}
    for s in sizes: hist[s] = hist.get(s, 0) + 1
    print("%s: V=%d E=%d F=%d genus=%s mesh_errors=%d A_viol=%d hull_faces=%s irregular_hull_faces=%d hull_edges=%d (double-counted) missing=%d" % (
        path, len(v), E, len(faces), (2 - chi) / 2, bad, aviol, dict(sorted(hist.items())), irregular, nedges, edges_missing))

for p in sys.argv[1:]:
    check(p)
