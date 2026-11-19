"""Query a window library written by the `windows` binary.

  python scripts/windows.py DIR list [--n N] [--min-n N] [--solid NAME] [--no-solid NAME] [--pure] [--rot K]
  python scripts/windows.py DIR show ID
  python scripts/windows.py DIR off ID OUT.off     example rim as an OFF file (window plane z = 0)

A window is a convex, planar, unit-edged polygon; `options[i]` lists the states (rows of
states.tsv: solid, rim edge, tilt) a piece on rim edge i can be in. --pure keeps windows whose
example rim has no solid face lying in the window plane.
"""
import sys, json, collections

def load(d):
    states = {}
    for l in open(d + '/states.tsv').read().split('\n')[1:]:
        if l:
            p = l.split('\t')
            states[int(p[0])] = dict(solid=p[1], a=int(p[2]), b=int(p[3]), tilt=float(p[4]), span=float(p[5]), xf=[float(x) for x in p[6].split()])
    wins = [json.loads(l) for l in open(d + '/windows.jsonl') if l.strip()]
    return states, wins

def solids_of(w, states):
    return [sorted({states[s]['solid'] for s in o}) for o in w['options']]

def turns_txt(w):
    c = collections.Counter(w['turns'])
    return ' '.join('%gx%.4f' % (k, t) for t, k in sorted(c.items()))

def main():
    d, cmd, rest = sys.argv[1], sys.argv[2], sys.argv[3:]
    states, wins = load(d)
    opt = lambda k: rest[rest.index(k) + 1] if k in rest else None
    opts = lambda k: [rest[i + 1] for i, x in enumerate(rest) if x == k]
    if cmd == 'list':
        print('id\tn\toutline\tsym\tinradius\trims\tverified\tfaces_in_plane\tturns\texample')
        for w in wins:
            every = {x for o in solids_of(w, states) for x in o}
            if opt('--n') and w['n'] != int(opt('--n')): continue
            if opt('--min-n') and w['n'] < int(opt('--min-n')): continue
            if opt('--rot') and w['rot'] != int(opt('--rot')): continue
            if any(s not in every for s in opts('--solid')): continue
            if any(s in every for s in opts('--no-solid')): continue
            ex = w['example']
            fip = ex['faces_in_plane'] if ex else -1
            if '--pure' in rest and fip != 0: continue
            sym = ('D' if w['mirror'] else 'C') + str(w['rot'])
            exs = ', '.join('%d %s' % (k, s) for s, k in sorted(collections.Counter(p['solid'] for p in ex['pieces']).items())) if ex else ''
            print('%d\t%d\t%d\t%s\t%.3f\t%.3g\t%s\t%d\t%s\t%s' % (w['id'], w['n'], w['outline'], sym, w['inradius'], w['rims'], w['verified'], fip, turns_txt(w), exs))
    elif cmd == 'show':
        w = wins[int(rest[0])]
        print('window %d: %d-gon, outline %d, rotations %d, mirror %s, rim walks %g' % (w['id'], w['n'], w['outline'], w['rot'], w['mirror'], w['rims']))
        ex = w['example']
        for i in range(w['n']):
            names = collections.Counter(states[s]['solid'] for s in w['options'][i])
            exs = states[ex['states'][i]] if ex else None
            print(' edge %2d  turn after %8.4f  example %-40s options %s' % (i, w['turns'][i], ('%s tilt %.3f' % (exs['solid'], exs['tilt'])) if exs else '-', dict(names)))
    elif cmd == 'off':
        import numpy as np
        w = wins[int(rest[0])]
        sol = {}
        L = open(opt('--solids') or 'solids.txt').read().split('\n'); i = 0
        while i < len(L):
            p = L[i].split()
            if not p: i += 1; continue
            nv, nf = int(p[2]), int(p[3])
            sol[p[1]] = (np.array([[float(x) for x in L[i + 1 + k].split()] for k in range(nv)]), [[int(x) for x in L[i + 1 + nv + k].split()[1:]] for k in range(nf)])
            i += 1 + nv + nf
        V, F = [], []
        for pc in w['example']['pieces']:
            v, f = sol[pc['solid']]
            x = np.array(pc['xf']); R, t = x[:9].reshape(3, 3), x[9:]
            base = len(V)
            V += list(v @ R.T + t)
            F += [[base + k for k in face] for face in f]
        with open(rest[1], 'w') as o:
            o.write('OFF\n%d %d 0\n' % (len(V), len(F)))
            for p in V: o.write('%.9f %.9f %.9f\n' % tuple(p))
            for f in F: o.write('%d %s\n' % (len(f), ' '.join(map(str, f))))
    else:
        sys.exit(__doc__)

main()
