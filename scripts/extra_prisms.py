"""Append k-gonal prisms and antiprisms to a copy of solids.txt (for window searches that want
face sizes the named solids do not have).  python scripts/extra_prisms.py solids.txt OUT 7 9 11 12"""
import math, sys

def prism(m):
    R = 0.5 / math.sin(math.pi / m)
    v = [(R * math.cos(2 * math.pi * i / m), R * math.sin(2 * math.pi * i / m), z) for z in (-0.5, 0.5) for i in range(m)]
    faces = [list(range(m))[::-1], list(range(m, 2 * m))]
    for i in range(m):
        j = (i + 1) % m
        faces.append([i, j, m + j, m + i])
    return v, faces

def antiprism(m):
    R = 0.5 / math.sin(math.pi / m)
    h = math.sqrt(1 - 0.25 / math.cos(math.pi / (2 * m)) ** 2)
    v = [(R * math.cos(2 * math.pi * i / m), R * math.sin(2 * math.pi * i / m), -h / 2) for i in range(m)]
    v += [(R * math.cos(2 * math.pi * (i + 0.5) / m), R * math.sin(2 * math.pi * (i + 0.5) / m), h / 2) for i in range(m)]
    faces = [list(range(m))[::-1], list(range(m, 2 * m))]
    for i in range(m):
        j = (i + 1) % m
        faces.append([i, j, m + i])
        faces.append([m + j, m + i, j])
    return v, faces

src, out, ks = sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3:]]
with open(out, 'w') as o:
    o.write(open(src).read().rstrip('\n') + '\n')
    for k in ks:
        for name, (v, f) in (('%d-prism' % k, prism(k)), ('%d-antiprism' % k, antiprism(k))):
            o.write('solid %s %d %d\n' % (name, len(v), len(f)))
            for p in v: o.write('%.17g %.17g %.17g\n' % p)
            for face in f: o.write('%d %s\n' % (len(face), ' '.join(map(str, face))))
