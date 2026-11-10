"""Build viewer/index.html. Usage: build_viewer.py picks.json
picks: [{"run": dir, "id": int, "name": str, "note": str, "grp": str}]"""
import csv, json, os, sys
import numpy as np
from scipy.spatial import ConvexHull

here = os.path.dirname(os.path.abspath(__file__))

def load_off(path):
    t = open(path).read().split()
    nv, nf = int(t[1]), int(t[2])
    v = np.array(t[4:4 + 3 * nv], dtype=float).reshape(nv, 3)
    faces, i = [], 4 + 3 * nv
    for _ in range(nf):
        n = int(t[i]); faces.append([int(x) for x in t[i + 1:i + 1 + n]]); i += n + 1
    return v, faces

def hull_faces(v):
    hull = ConvexHull(v)
    planes = {}
    for s, eq in zip(hull.simplices, hull.equations):
        planes.setdefault(tuple(np.round(eq, 5)), []).append(s)
    out = []
    for key, tris in planes.items():
        n = np.array(key[:3])
        ids = sorted({int(i) for s in tris for i in s})
        p = v[ids]; c = p.mean(axis=0)
        e1 = p[0] - c; e1 /= np.linalg.norm(e1); e2 = np.cross(n, e1)
        ang = np.arctan2((p - c) @ e2, (p - c) @ e1)
        order = [ids[i] for i in np.argsort(ang)]
        m = len(order)
        corners = [order[i] for i in range(m) if np.linalg.norm(np.cross(v[order[i]] - v[order[i - 1]], v[order[(i + 1) % m]] - v[order[i]])) > 1e-7]
        pc = v[corners]; r = np.linalg.norm(pc - pc.mean(axis=0), axis=1)
        out.append(dict(i=corners, r=int(r.max() - r.min() < 1e-6)))
    return out

picks = json.load(open(sys.argv[1]))
items = []
for p in picks:
    row = next(r for r in csv.DictReader(open(os.path.join(p["run"], "results.tsv")), delimiter="\t") if int(r["id"]) == p["id"])
    off = os.path.join(p["run"], "off", "%06d.off" % p["id"])
    v, faces = load_off(off)
    v = v - v.mean(axis=0)
    names = open(off + ".pieces").read().split("\n")
    content = [c.split("x", 1) for c in row["content"].split("+")]
    pieces = [c[1] for c in content]
    hf = hull_faces(v)
    hist = {}
    for h in hf:
        k = "%d%s" % (len(h["i"]), "" if h["r"] else "*")
        hist[k] = hist.get(k, 0) + 1
    R = float(np.linalg.norm(v, axis=1).max())
    items.append(dict(name=p["name"], note=p["note"], grp=p.get("grp", ""), genus=int(row["genus"]),
        sym={"60": "icosahedral", "24": "octahedral", "12": "tetrahedral"}[row["group"]], npieces=int(row["pieces"]), F=len(faces),
        hull=", ".join("%d × %s-gon" % (c, k.replace("*", " irregular")) if False else "%d×%s" % (c, k) for k, c in sorted(hist.items(), key=lambda x: int(x[0].rstrip("*")))),
        round="%.3f" % float(row["roundness"]), diam="%.1f" % (2 * R), cls="(R)(A)(Q)(T)" if row["class"] == "RAQT" else "near miss, %s hull edges missing" % row["Q_missing"],
        chain=row["chain"], R=round(R, 3), pieces=[n.replace("-", " ") for n in pieces], counts=[int(c[0]) for c in content],
        v=[round(float(x), 4) for x in v.reshape(-1)], f=faces, p=[pieces.index(n) for n in names], h=hf))
    print(p["name"], row["content"], hist)
html = open(os.path.join(here, "template.html")).read().replace("/*DATA*/[]", json.dumps(items, separators=(",", ":")))
open(os.path.join(here, "index.html"), "w").write(html)
print("wrote index.html", len(html) // 1024, "KB")
