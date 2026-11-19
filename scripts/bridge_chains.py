"""Which chains of rim pieces can bridge the two antiprism ends of Webb's 12-gon window?

  python scripts/bridge_chains.py FULL_DIR

FULL_DIR is a window run made with --full-tables (all usable states and corners, about 30 MB),
e.g. `windows --min-n 20 --full-tables --out /tmp/full`. A C2 window needs the chain's turns,
both end corners included, to add up to 126.87 degrees.
"""
import json, sys, collections, math
full = sys.argv[1]
st = {}
for l in open(full + '/states.tsv').read().split('\n')[1:]:
    if l:
        p = l.split('\t'); st[int(p[0])] = (p[1], float(p[4]), float(p[5]))
out = collections.defaultdict(list)
for l in open(full + '/corners.tsv').read().split('\n')[1:]:
    if l:
        p = l.split('\t'); out[int(p[0])].append((int(p[1]), float(p[2]), p[5]))
for k in out: out[k] = sorted(set(out[k]))
w = [json.loads(l) for l in open('windows/open/windows.jsonl')][27]
ex = w['example']['states']
print('webb rim:', [(st[s][0], round(st[s][1], 3)) for s in ex])
s_in, rot, s_out = ex[0], ex[1], ex[2]
assert st[s_in][0] == 'pentagonal-antiprism' and st[rot][0] == 'pentagonal-rotunda'
print('S_out state', s_out, st[s_out], 'S_in', s_in, st[s_in])
print('corners leaving S_out:')
for to, t, via in out[s_out]: print('   turn %.4f -> %s tilt %.3f/%.3f via %s' % (t, st[to][0], st[to][1], st[to][2], via))
inn = [(a, t, via) for a in out for (to, t, via) in out[a] if to == s_in]
print('corners entering S_in:')
for a, t, via in inn: print('   turn %.4f <- %s tilt %.3f/%.3f via %s' % (t, st[a][0], st[a][1], st[a][2], via))
need = 180 - 2 * 26.565051
print('chain turn sum needed for a C2 window: %.4f' % need)
for k in (3, 4, 5, 6):
    res = []
    cnt = 0
    def go(cur, path, tot):
        global cnt
        if len(path) == k:
            for to, t, via in out[cur]:
                if to == s_in:
                    cnt += 1; res.append((tot + t, path[:]))
            return
        for to, t, via in out[cur]:
            if tot + t > 2 * need: continue
            path.append((to, t)); go(to, path, tot + t); path.pop()
    go(s_out, [], 0.0)
    sums = collections.Counter(round(r[0], 3) for r in res)
    print('k = %d pieces between the antiprisms: %d chains, %d distinct turn sums' % (k, cnt, len(sums)))
    near = sorted(sums.items(), key=lambda x: abs(x[0] - need))[:6]
    print('   closest sums to %.3f:' % need, near)
    exact = [r for r in res if abs(r[0] - need) < 1e-3]
    names = collections.Counter(tuple(st[s][0] for s, t in r[1]) for r in exact)
    for nm, c in names.most_common(8): print('   exact:', c, nm)
    # unequal pairs: sums a + b = 2 * need
    pairs = [(a, b) for a in sums for b in sums if a < b and abs(a + b - 2 * need) < 1e-3]
    print('   unequal pairs of sums adding to %.3f: %s' % (2 * need, pairs[:10]))
