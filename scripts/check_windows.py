"""Independent check of a window library (numpy + scipy, shares no code with the Rust search).

  python scripts/check_windows.py DIR [--solids solids.txt]

For the example rim of every window: rebuilds the pieces from solids.txt, then checks that
nothing lies above the window plane, that the points in the plane have the stored polygon as
their convex hull, that the polygon is convex with unit edges, that every polygon edge is an edge
of a piece, and that no two pieces overlap (largest ball inside the intersection, by LP).
"""
import sys, json, itertools
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull

d = sys.argv[1]
sp = sys.argv[sys.argv.index('--solids') + 1] if '--solids' in sys.argv else 'solids.txt'
sol = {}
L = open(sp).read().split('\n'); i = 0
while i < len(L):
    p = L[i].split()
    if not p: i += 1; continue
    nv, nf = int(p[2]), int(p[3])
    sol[p[1]] = (np.array([[float(x) for x in L[i + 1 + k].split()] for k in range(nv)]), [[int(x) for x in L[i + 1 + nv + k].split()[1:]] for k in range(nf)])
    i += 1 + nv + nf

def halfspaces(v):
    h = ConvexHull(v)
    return np.unique(np.round(h.equations, 9), axis=0)

def overlap(ha, hb):
    """radius of the largest ball inside both solids (<= 0: interiors disjoint)"""
    A = np.vstack([ha, hb])
    c = np.zeros(4); c[3] = -1
    r = linprog(c, A_ub=np.hstack([A[:, :3], np.ones((len(A), 1))]), b_ub=-A[:, 3], bounds=[(None, None)] * 3 + [(None, 1)])
    return r.x[3] if r.status == 0 else -1

bad = 0; n_ok = 0; worst = 0.0
for line in open(d + '/windows.jsonl'):
    w = json.loads(line)
    ex = w['example']
    if not ex: continue
    errs = []
    poly = np.array(w['poly']); n = len(poly)
    pcs = []
    for pc in ex['pieces']:
        v, f = sol[pc['solid']]
        x = np.array(pc['xf']); R, t = x[:9].reshape(3, 3), x[9:]
        if abs(np.linalg.det(R) - 1) > 1e-9: errs.append('not a rotation')
        pcs.append((v @ R.T + t, f))
    allv = np.vstack([v for v, f in pcs])
    if allv[:, 2].max() > 1e-7: errs.append('above plane %.2e' % allv[:, 2].max())
    flat = allv[allv[:, 2] > -1e-7][:, :2]
    flat = np.unique(np.round(flat, 7), axis=0)  # coincident corners of neighbouring pieces
    hv = flat[ConvexHull(flat).vertices]
    # hull corners must be exactly the polygon corners
    if len(hv) != n or max(np.abs(poly - q).sum(1).min() for q in hv) > 1e-6: errs.append('hull of in-plane points is not the polygon (%d vs %d)' % (len(hv), n))
    side = np.linalg.norm(np.roll(poly, -1, 0) - poly, axis=1)
    if np.abs(side - 1).max() > 1e-7: errs.append('edge length')
    e = np.roll(poly, -1, 0) - poly; e2 = np.roll(e, -1, 0)
    if (e[:, 0] * e2[:, 1] - e[:, 1] * e2[:, 0]).min() < 1e-6: errs.append('not strictly convex')
    # each polygon edge is an edge of some piece
    for k in range(n):
        a = np.append(poly[k], 0); b = np.append(poly[(k + 1) % n], 0)
        hit = False
        for v, f in pcs:
            for face in f:
                for m in range(len(face)):
                    p, q = v[face[m]], v[face[(m + 1) % len(face)]]
                    if np.abs(p - a).max() < 1e-6 and np.abs(q - b).max() < 1e-6: hit = True
                    if np.abs(p - b).max() < 1e-6 and np.abs(q - a).max() < 1e-6: hit = True
        if not hit: errs.append('polygon edge %d is not a piece edge' % k)
    hs = [halfspaces(v) for v, f in pcs]
    for i, j in itertools.combinations(range(len(pcs)), 2):
        r = overlap(hs[i], hs[j])
        worst = max(worst, r)
        if r > 1e-6: errs.append('pieces %d and %d overlap (ball %.2e)' % (i, j, r))
    if errs: bad += 1; print('window', w['id'], 'n', n, errs)
    else: n_ok += 1
print('checked %d windows: %d pass, %d fail; deepest overlap ball %.2e' % (n_ok + bad, n_ok, bad, worst))
