"""strip.py RUN ID : rebuild all pieces of a result, drop prisms/antiprisms, compare hulls."""
import sys, json, collections
import numpy as np
from scipy.spatial import ConvexHull
run, rid = sys.argv[1], int(sys.argv[2])
solids = {}
L = open('solids.txt').read().split('\n'); i = 0
while i < len(L):
    p = L[i].split()
    if not p: i += 1; continue
    nv, nf = int(p[2]), int(p[3])
    v = np.array([[float(x) for x in L[i+1+k].split()] for k in range(nv)])
    f = [[int(x) for x in L[i+1+nv+k].split()[1:]] for k in range(nf)]
    solids[p[1]] = (v, f); i += 1 + nv + nf
rec = next(r for r in map(json.loads, open(run + '/results.jsonl')) if r['id'] == rid)
def xf(a): a = np.array(a); return a[:9].reshape(3, 3), a[9:]
gens = [xf(g) for g in rec['gens']]
grp = [(np.eye(3), np.zeros(3))]
k = 0
while k < len(grp):
    for (R, t) in gens:
        n = (R @ grp[k][0], R @ grp[k][1] + t)
        if not any(np.abs(n[0] - e[0]).max() < 1e-6 for e in grp): grp.append(n)
    k += 1
centre = np.array(rec['centre'])
pieces = []  # (name, verts)
for rp in rec['reps']:
    R0, t0 = xf(rp['xf']); v, f = solids[rp['solid']]
    for (R, t) in grp:
        c = R @ t0 + t
        if any(n == rp['solid'] and np.linalg.norm(c - cc) < 1e-5 for n, cc, _ in pieces): continue
        pieces.append((rp['solid'], c, (v @ R0.T + t0) @ R.T + t - centre))
print('group', len(grp), 'pieces', dict(collections.Counter(n for n, _, _ in pieces)))
strut = lambda n: n.endswith('-prism') or n.endswith('-antiprism')
def hull_faces(P):
    h = ConvexHull(P); planes = {}
    for s, eq in zip(h.simplices, h.equations): planes.setdefault(tuple(np.round(eq, 5)), set()).update(int(x) for x in s)
    out = []
    for key, ids in planes.items():
        n = np.array(key[:3]); ids = sorted(ids); p = P[ids]; c = p.mean(0)
        e1 = p[0] - c; e1 /= np.linalg.norm(e1); e2 = np.cross(n, e1)
        o = [ids[j] for j in np.argsort(np.arctan2((p - c) @ e2, (p - c) @ e1))]
        m = len(o)
        cor = [o[j] for j in range(m) if np.linalg.norm(np.cross(P[o[j]] - P[o[j-1]], P[o[(j+1) % m]] - P[o[j]])) > 1e-7]
        sides = [np.linalg.norm(P[cor[j]] - P[cor[j-1]]) for j in range(len(cor))]
        out.append((len(cor), tuple(sorted(round(s, 4) for s in sides)), -key[3]))
    return out
def report(label, sel):
    P = np.vstack([v for n, c, v in pieces if sel(n)])
    P = np.unique(np.round(P, 6), axis=0)
    hf = hull_faces(P)
    print(f'-- {label}: {len(P)} distinct vertices, circumradius {np.linalg.norm(P, axis=1).max():.4f}')
    for (n, sides, d), cnt in sorted(collections.Counter(hf).items()):
        print(f'   {cnt:3d} x {n}-gon at distance {d:.4f} sides {dict(collections.Counter(sides))}')
    return P
full = report('full toroid', lambda n: True)
nodes = report('nodes only (prisms and antiprisms removed)', lambda n: not strut(n))
hv = lambda P: {tuple(x) for x in np.round(P[ConvexHull(P).vertices], 5)}
a, b = hv(full), hv(nodes)
print('hull vertices: full', len(a), 'nodes-only', len(b), 'shared', len(a & b))
for n, c, v in pieces[:0]: pass
rad = collections.defaultdict(set)
for n, c, v in pieces: rad[n].add(round(float(np.linalg.norm(c - centre)), 4))
print('piece centre distances', {k: sorted(v) for k, v in rad.items()})
