//! Window library: convex, planar, unit-edged n-gons that can be an open face of the convex hull,
//! together with the ring of solids whose edges form the rim. See windows/README.md.
//!
//! Model. A rim piece is a solid lying under the window plane and touching it along one edge.
//! Its state is (solid, directed edge up to the solid's rotations, tilt of the plane about the edge).
//! A corner is two rim pieces glued face to face whose rim edges meet at a vertex of the shared
//! face and span a plane that supports both. Corners are arcs between states; a window is a closed
//! walk whose turns add up to 360 degrees and whose unit steps close up.
#![allow(dead_code)]
#[path = "../geom.rs"]
mod geom;
#[path = "../simd.rs"]
mod simd;
#[path = "../solids.rs"]
mod solids;

use geom::*;
use solids::Solid;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::io::Write;

const TOL: f64 = 1e-9;
const TILT_BIN: f64 = 1e-6;

struct Sym { vperm: Vec<Vec<usize>> }

struct State { solid: usize, a: usize, b: usize, tilt: f64, span: f64, xf: Xf }

#[derive(Clone)]
struct Arc { from: u32, to: u32, turn: f64, face: usize, aplanar: bool,
    /// hidden piece joining the two rim pieces at the corner: solid and placement in the rim frame of `from`
    via: Option<(usize, Xf)> }

struct Graph { states: Vec<State>, arcs: Vec<Arc>, index: HashMap<(usize, usize, usize, i64), u32> }

fn vperms(s: &Solid) -> Vec<Vec<usize>> {
    s.syms.iter().map(|r| s.v.iter().map(|p| s.v.iter().position(|q| dist(*q, mv(r, *p)) < 1e-6).unwrap()).collect()).collect()
}

/// the face that has a -> b in its counter-clockwise boundary
fn left_face(s: &Solid, a: usize, b: usize) -> usize {
    for (fi, f) in s.faces.iter().enumerate() {
        for i in 0..f.len() { if f[i] == a && f[(i + 1) % f.len()] == b { return fi; } }
    }
    panic!("no face");
}

/// rim frame: a at the origin, b at (1,0,0), plane z = 0, solid below
fn rim_xf(s: &Solid, a: usize, b: usize, n: V3) -> Xf {
    let x = unit(sub(s.v[b], s.v[a]));
    let y = cross(n, x);
    let r = [x, y, n];
    Xf { r, t: scale(mv(&r, s.v[a]), -1.0) }
}

impl Graph {
    /// state of solid si touching the plane with unit normal n (solid coordinates, pointing away
    /// from the solid) along the directed edge a -> b
    fn state(&mut self, lib: &[Solid], sy: &[Sym], si: usize, a: usize, b: usize, n: V3) -> u32 {
        let s = &lib[si];
        let nl = s.fnorm[left_face(s, a, b)];
        let nr = s.fnorm[left_face(s, b, a)];
        let tilt = dot(n, nl).clamp(-1.0, 1.0).acos();
        let span = dot(nr, nl).clamp(-1.0, 1.0).acos();
        let (ra, rb) = sy[si].vperm.iter().map(|p| (p[a], p[b])).min().unwrap();
        let k = (tilt / TILT_BIN).round() as i64;
        for d in [0, -1, 1] {
            if let Some(&id) = self.index.get(&(si, ra, rb, k + d)) { return id; }
        }
        let xf = rim_xf(s, a, b, n);
        let id = self.states.len() as u32;
        self.states.push(State { solid: si, a, b, tilt, span, xf });
        self.index.insert((si, ra, rb, k), id);
        id
    }
}

fn build(lib: &[Solid], active: &[bool], via: bool) -> Graph {
    let sy: Vec<Sym> = lib.iter().map(|s| Sym { vperm: vperms(s) }).collect();
    let nbr: Vec<Vec<Vec<usize>>> = lib.iter().map(|s| {
        let mut n = vec![Vec::new(); s.v.len()];
        for &(a, b) in &s.edges { n[a].push(b); n[b].push(a); }
        n
    }).collect();
    let mut g = Graph { states: Vec::new(), arcs: Vec::new(), index: HashMap::new() };
    let mut seen: HashMap<(u32, u32, i64), ()> = HashMap::new();
    let flip = |p: V3| -> V3 { [p[0], -p[1], -p[2]] };
    let (mut straight, mut tested) = (0u64, 0u64);
    for ai in 0..lib.len() {
        if !active[ai] { continue; }
        let sa = &lib[ai];
        for &(fa, ka) in &sa.flags {
            let m = sa.faces[fa].len();
            let xa = sa.frames[fa][ka].inv();
            let pa: Vec<V3> = sa.v.iter().map(|p| xa.ap(*p)).collect();
            let va = sa.faces[fa][ka];
            let o = pa[va];
            // normals of A's faces next to the glue face, per glue edge, for the (A) test
            for bi in ai..lib.len() {
                if !active[bi] { continue; }
                let sb = &lib[bi];
                for &(fb, kb) in &sb.flags {
                    if sb.faces[fb].len() != m { continue; }
                    if bi == ai && (fb, kb) < (fa, ka) { continue; }
                    let xb = sb.frames[fb][kb].inv();
                    let pb: Vec<V3> = sb.v.iter().map(|p| flip(xb.ap(*p))).collect();
                    let vb = sb.faces[fb][kb];
                    // (A): no face of A coplanar with the face of B across a glue edge
                    let mut aplanar = true;
                    for i in 0..m {
                        let (a0, a1) = (sa.faces[fa][(ka + i) % m], sa.faces[fa][(ka + i + 1) % m]);
                        let na = mv(&xa.r, sa.fnorm[left_face(sa, a1, a0)]);
                        // B's face is mirrored: vertex kb - i sits on A's vertex ka + i
                        let (b0, b1) = (sb.faces[fb][(kb + m - i) % m], sb.faces[fb][(kb + m - i - 1) % m]);
                        let nb = flip(mv(&xb.r, sb.fnorm[left_face(sb, b0, b1)]));
                        debug_assert!(dist(pa[a0], pb[b0]) < 1e-6 && dist(pa[a1], pb[b1]) < 1e-6);
                        if dot(na, nb) > 1.0 - 1e-9 { aplanar = false; }
                    }
                    for &ea in &nbr[ai][va] {
                        let da = sub(pa[ea], o);
                        for &eb in &nbr[bi][vb] {
                            let db = sub(pb[eb], o);
                            let c = cross(da, db);
                            let cn = norm(c);
                            if cn < 1e-7 { if dot(da, db) < 0.0 { straight += 1; } continue; }
                            tested += 1;
                            let mut n = scale(c, 1.0 / cn);
                            // quick sign and support test on the neighbours of the corner vertex
                            let mut sign = 0i32;
                            let mut ok = true;
                            for p in nbr[ai][va].iter().map(|&i| pa[i]).chain(nbr[bi][vb].iter().map(|&i| pb[i])) {
                                let h = dot(n, sub(p, o));
                                if h > TOL { if sign < 0 { ok = false; break; } sign = 1; }
                                if h < -TOL { if sign > 0 { ok = false; break; } sign = -1; }
                            }
                            if !ok || sign == 0 { continue; }
                            if sign > 0 { n = scale(n, -1.0); }
                            if pa.iter().chain(pb.iter()).any(|p| dot(n, sub(*p, o)) > TOL) { continue; }
                            // interior angle between the two rim edges; walking counter-clockwise seen from +n
                            let theta = dot(da, db).clamp(-1.0, 1.0).acos();
                            let turn = PI - theta;
                            let na = mv(&mt(&xa.r), n);
                            let nb = mv(&mt(&xb.r), flip(n));
                            let a_first = dot(cross(db, da), n) > 0.0;
                            let (from, to) = if a_first {
                                (g.state(lib, &sy, ai, ea, va, na), g.state(lib, &sy, bi, vb, eb, nb))
                            } else {
                                (g.state(lib, &sy, bi, eb, vb, nb), g.state(lib, &sy, ai, va, ea, na))
                            };
                            let key = (from, to, (turn / TILT_BIN).round() as i64);
                            if seen.insert(key, ()).is_none() {
                                g.arcs.push(Arc { from, to, turn, face: m, aplanar, via: None });
                            } else if aplanar {
                                // the same corner reached through another gluing: keep the better flag
                                if let Some(x) = g.arcs.iter_mut().rev().find(|x| x.via.is_none() && x.from == from && x.to == to && (x.turn - turn).abs() < 2.0 * TILT_BIN) { x.aplanar = true; }
                            }
                        }
                    }
                }
            }
        }
    }
    eprintln!("corner tests {} straight (180 degree) pairs skipped {}", tested, straight);
    if via { build_via(lib, active, &sy, &nbr, &mut g); }
    // a piece with a face in the window plane: the rim can follow that face round to its next edge
    let mut i = 0;
    while i < g.states.len() {
        let (si, a, b, tilt) = (g.states[i].solid, g.states[i].a, g.states[i].b, g.states[i].tilt);
        if tilt < 1e-7 {
            let s = &lib[si];
            let fi = left_face(s, a, b);
            let f = &s.faces[fi];
            let k = f.iter().position(|&x| x == b).unwrap();
            let c = f[(k + 1) % f.len()];
            let to = g.state(lib, &sy, si, b, c, s.fnorm[fi]);
            g.arcs.push(Arc { from: i as u32, to, turn: 2.0 * PI / f.len() as f64, face: 0, aplanar: true, via: None });
        }
        i += 1;
    }
    g
}

const TWO_PI: f64 = 2.0 * PI;
const CLOSE: f64 = 1e-7;

struct Piece { solid: usize, state: u32, r: M3, v: Vec<V3>, c: V3 }

/// largest separation over the usual axes; negative means the solids overlap by that much
fn sat(lib: &[Solid], a: &Piece, b: &Piece) -> f64 {
    let gap = |ax: V3| -> f64 {
        let (mut a0, mut a1, mut b0, mut b1) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for p in &a.v { let h = dot(ax, *p); a0 = a0.min(h); a1 = a1.max(h); }
        for p in &b.v { let h = dot(ax, *p); b0 = b0.min(h); b1 = b1.max(h); }
        (b0 - a1).max(a0 - b1)
    };
    let mut best = f64::MIN;
    for n in &lib[a.solid].fnorm { best = best.max(gap(mv(&a.r, *n))); if best > -1e-7 { return best; } }
    for n in &lib[b.solid].fnorm { best = best.max(gap(mv(&b.r, *n))); if best > -1e-7 { return best; } }
    for ea in &lib[a.solid].edirs {
        let da = mv(&a.r, *ea);
        for eb in &lib[b.solid].edirs {
            let c = cross(da, mv(&b.r, *eb));
            if norm(c) < 1e-9 { continue; }
            best = best.max(gap(unit(c)));
            if best > -1e-7 { return best; }
        }
    }
    best
}

/// One side of a corner: a solid glued by face f (vertex k on the corner) with rim edge va -> e.
struct Half { solid: usize, va: usize, e: usize, m: usize, xa: Xf, pts: Vec<V3>, d: V3, nb: Vec<V3>, n0: V3, w: V3, span: f64 }

fn fm(p: V3) -> V3 { [p[0], -p[1], -p[2]] }
fn fm_xf(x: &Xf) -> Xf { Xf { r: [x.r[0], scale(x.r[1], -1.0), scale(x.r[2], -1.0)], t: fm(x.t) } }

/// can a plane through the half's rim edge support the half and stay off the directions q?
fn half_fits(h: &Half, q: &[V3]) -> bool {
    let (mut lo, mut hi) = (0.0f64, h.span);
    for qj in q {
        let (a, b) = (dot(h.n0, *qj), dot(h.w, *qj));
        let g1 = a * hi.cos() + b * hi.sin();
        let g0 = a * lo.cos() + b * lo.sin();
        if g0 <= TOL && g1 <= TOL { continue; }
        if g0 > TOL && g1 > TOL { return false; }
        let mut t = (-a).atan2(b);
        while t < lo - 1e-12 { t += PI; }
        while t > hi + 1e-12 { t -= PI; }
        if t < lo - 1e-9 { return false; }
        if g0 > TOL { lo = t; } else { hi = t; }
        if lo > hi + 1e-9 { return false; }
    }
    true
}

/// Corners where the two rim pieces are both glued to a third piece that touches the window plane
/// only at the corner vertex.
fn build_via(lib: &[Solid], active: &[bool], sy: &[Sym], nbr: &[Vec<Vec<usize>>], g: &mut Graph) {
    let mut halves: Vec<Half> = Vec::new();
    let mut by_m: HashMap<usize, Vec<usize>> = HashMap::new();
    for ai in 0..lib.len() {
        if !active[ai] { continue; }
        let sa = &lib[ai];
        for &(fa, ka) in &sa.flags {
            let xa = sa.frames[fa][ka].inv();
            let pts: Vec<V3> = sa.v.iter().map(|p| xa.ap(*p)).collect();
            let va = sa.faces[fa][ka];
            for &e in &nbr[ai][va] {
                let n0 = mv(&xa.r, sa.fnorm[left_face(sa, va, e)]);
                let n1 = mv(&xa.r, sa.fnorm[left_face(sa, e, va)]);
                let span = dot(n0, n1).clamp(-1.0, 1.0).acos();
                let w = unit(sub(n1, scale(n0, dot(n0, n1))));
                by_m.entry(sa.faces[fa].len()).or_default().push(halves.len());
                halves.push(Half { solid: ai, va, e, m: sa.faces[fa].len(), xa, pts: pts.clone(), d: sub(pts[e], pts[va]), nb: nbr[ai][va].iter().map(|&i| sub(pts[i], pts[va])).collect(), n0, w, span });
            }
        }
    }
    let mut seen: std::collections::HashSet<(u32, u32, i64, usize, [i64; 3])> = Default::default();
    let (mut tested, mut found) = (0u64, 0u64);
    let empty = Vec::new();
    for ci in 0..lib.len() {
        if !active[ci] { continue; }
        let sc = &lib[ci];
        for &(f1, k1) in &sc.flags {
            let vc = sc.faces[f1][k1];
            let t1 = fm_xf(&sc.frames[f1][k1].inv());
            let o = t1.ap(sc.v[vc]);
            let cq: Vec<V3> = nbr[ci][vc].iter().map(|&i| sub(t1.ap(sc.v[i]), o)).collect();
            let cpts: Vec<V3> = sc.v.iter().map(|p| t1.ap(*p)).collect();
            let ha: Vec<usize> = by_m.get(&sc.faces[f1].len()).unwrap_or(&empty).iter().cloned().filter(|&h| half_fits(&halves[h], &cq)).collect();
            if ha.is_empty() { continue; }
            for f2 in 0..sc.faces.len() {
                if f2 == f1 { continue; }
                let k2 = match sc.faces[f2].iter().position(|&x| x == vc) { Some(k) => k, None => continue };
                let t2 = fm_xf(&sc.frames[f2][k2].inv());
                let w2 = t1.mul(&t2.inv());
                let w2r_t = mt(&w2.r);
                let cq2: Vec<V3> = cq.iter().map(|q| mv(&w2r_t, *q)).collect();
                for &hb in by_m.get(&sc.faces[f2].len()).unwrap_or(&empty) {
                    let b = &halves[hb];
                    if !half_fits(b, &cq2) { continue; }
                    let db = mv(&w2.r, b.d);
                    let bnb: Vec<V3> = b.nb.iter().map(|q| mv(&w2.r, *q)).collect();
                    let mut bpts: Option<Vec<V3>> = None;
                    for &hai in &ha {
                        let a = &halves[hai];
                        let c = cross(a.d, db);
                        let cn = norm(c);
                        if cn < 1e-7 { continue; }
                        tested += 1;
                        let mut n = scale(c, 1.0 / cn);
                        let mut sign = 0i32;
                        let mut ok = true;
                        for q in cq.iter().chain(a.nb.iter()).chain(bnb.iter()) {
                            let h = dot(n, *q);
                            if h > TOL { if sign < 0 { ok = false; break; } sign = 1; }
                            if h < -TOL { if sign > 0 { ok = false; break; } sign = -1; }
                        }
                        if !ok || sign == 0 { continue; }
                        if sign > 0 { n = scale(n, -1.0); }
                        let bp = bpts.get_or_insert_with(|| b.pts.iter().map(|p| w2.ap(*p)).collect());
                        if a.pts.iter().chain(cpts.iter()).chain(bp.iter()).any(|p| dot(n, sub(*p, o)) > TOL) { continue; }
                        // the hidden piece may touch the plane at the corner vertex only
                        if cpts.iter().filter(|p| dot(n, sub(**p, o)) > -1e-7).count() != 1 { continue; }
                        // the two rim pieces must not run into each other
                        let pa = Piece { solid: a.solid, state: 0, r: a.xa.r, v: a.pts.clone(), c: [0.0; 3] };
                        let xbw = w2.mul(&b.xa);
                        let pb = Piece { solid: b.solid, state: 0, r: xbw.r, v: bp.clone(), c: [0.0; 3] };
                        if sat(lib, &pa, &pb) < -1e-6 { continue; }
                        let theta = (dot(a.d, db) / (norm(a.d) * norm(db))).clamp(-1.0, 1.0).acos();
                        let turn = PI - theta;
                        let na = mv(&mt(&a.xa.r), n);
                        let nb = mv(&mt(&xbw.r), n);
                        let a_first = dot(cross(db, a.d), n) > 0.0;
                        let (from, to, fx) = if a_first {
                            (g.state(lib, sy, a.solid, a.e, a.va, na), g.state(lib, sy, b.solid, b.va, b.e, nb), rim_xf(&lib[a.solid], a.e, a.va, na).mul(&a.xa.inv()))
                        } else {
                            (g.state(lib, sy, b.solid, b.e, b.va, nb), g.state(lib, sy, a.solid, a.va, a.e, na), rim_xf(&lib[b.solid], b.e, b.va, nb).mul(&xbw.inv()))
                        };
                        let vx = fx.mul(&t1);
                        let mut cc = [0.0; 3];
                        for p in &sc.v { cc = add(cc, vx.ap(*p)); }
                        let key = (from, to, (turn / TILT_BIN).round() as i64, ci, [(cc[0] * 1e4 / sc.v.len() as f64).round() as i64, (cc[1] * 1e4 / sc.v.len() as f64).round() as i64, (cc[2] * 1e4 / sc.v.len() as f64).round() as i64]);
                        if seen.insert(key) {
                            found += 1;
                            g.arcs.push(Arc { from, to, turn, face: sc.faces[f1].len() * 100 + sc.faces[f2].len(), aplanar: true, via: Some((ci, vx)) });
                        }
                    }
                }
            }
        }
    }
    eprintln!("hidden-piece corner tests {} corners {}", tested, found);
}

/// The corner graph restricted to usable states, with turns grouped into classes.
struct Net {
    ns: usize,
    cls: Vec<f64>,
    /// out[state]: (turn class, next state), sorted
    out: Vec<Vec<(u16, u32)>>,
    cidx: Vec<u32>,
    na: usize,
    /// d[cidx[x] * na + cidx[s]]: smallest total turn of a walk x -> s
    d: Vec<f32>,
    /// graph arcs behind each (from, to, turn class), direct corners first
    between: HashMap<(u32, u32, u16), Vec<u32>>,
}

impl Net {
    fn step(&self, cur: u32, c: u16) -> impl Iterator<Item = u32> + '_ {
        let o = &self.out[cur as usize];
        let lo = o.partition_point(|x| x.0 < c);
        o[lo..].iter().take_while(move |x| x.0 == c).map(|x| x.1)
    }
    fn dist(&self, x: u32, s: u32) -> f64 { self.d[self.cidx[x as usize] as usize * self.na + self.cidx[s as usize] as usize] as f64 }
}

fn make_net(g: &Graph, coplanar: bool, no_faces: bool) -> Net {
    let ns = g.states.len();
    // a state whose right-hand face lies in the plane puts that face outside the polygon
    let mut alive: Vec<bool> = g.states.iter().map(|s| !(s.span - s.tilt < 1e-7) && !(no_faces && s.tilt < 1e-7)).collect();
    let ok = |a: &Arc, alive: &[bool]| alive[a.from as usize] && alive[a.to as usize] && (coplanar || a.aplanar);
    loop {
        let (mut i, mut o) = (vec![false; ns], vec![false; ns]);
        for a in &g.arcs { if ok(a, &alive) { o[a.from as usize] = true; i[a.to as usize] = true; } }
        let mut ch = false;
        for k in 0..ns { if alive[k] && !(i[k] && o[k]) { alive[k] = false; ch = true; } }
        if !ch { break; }
    }
    let mut ts: Vec<f64> = g.arcs.iter().filter(|a| ok(a, &alive)).map(|a| a.turn).collect();
    ts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut cls: Vec<f64> = Vec::new();
    for t in ts { if cls.last().map_or(true, |l| t - l > 1e-7) { cls.push(t); } }
    let class = |t: f64| -> u16 { let k = cls.partition_point(|c| *c < t - 1e-7); k as u16 };
    let mut out: Vec<Vec<(u16, u32)>> = vec![Vec::new(); ns];
    let mut radj: Vec<Vec<(u32, f64)>> = vec![Vec::new(); ns];
    let mut between: HashMap<(u32, u32, u16), Vec<u32>> = HashMap::new();
    for (k, a) in g.arcs.iter().enumerate() {
        if ok(a, &alive) {
            out[a.from as usize].push((class(a.turn), a.to));
            radj[a.to as usize].push((a.from, a.turn));
            let e = between.entry((a.from, a.to, class(a.turn))).or_default();
            if a.via.is_none() { e.insert(0, k as u32); } else { e.push(k as u32); }
        }
    }
    for o in out.iter_mut() { o.sort(); o.dedup(); }
    let mut cidx = vec![u32::MAX; ns];
    let mut na = 0;
    for k in 0..ns { if alive[k] { cidx[k] = na as u32; na += 1; } }
    let mut d = vec![f32::INFINITY; na * na];
    let mut heap = std::collections::BinaryHeap::new();
    let mut dd = vec![f64::INFINITY; ns];
    for s in 0..ns {
        if !alive[s] { continue; }
        for x in dd.iter_mut() { *x = f64::INFINITY; }
        dd[s] = 0.0;
        heap.push((std::cmp::Reverse(0u64), s as u32));
        while let Some((std::cmp::Reverse(k), u)) = heap.pop() {
            if k as f64 * 1e-12 > dd[u as usize] + 1e-11 { continue; }
            for &(v, t) in &radj[u as usize] {
                let nd = dd[u as usize] + t;
                if nd < dd[v as usize] - 1e-11 { dd[v as usize] = nd; heap.push((std::cmp::Reverse((nd * 1e12) as u64), v)); }
            }
        }
        // the walk back to s itself must be a real closed walk, not the empty one
        let mut selfd = f64::INFINITY;
        for &(v, t) in &radj[s] { selfd = selfd.min(dd[v as usize] + t); }
        let ci = cidx[s] as usize;
        for x in 0..ns { if alive[x] { d[cidx[x] as usize * na + ci] = (if x == s { selfd } else { dd[x] } - 1e-5) as f32; } }
    }
    Net { ns, cls, out, cidx, na, d, between }
}

struct ShapeSearch<'a> { net: &'a Net, nmin: usize, nmax: usize, t1: f64, c1: usize, seq: Vec<u16>, found: Vec<Vec<u16>>, nodes: u64, node_cap: u64, second: usize }

impl<'a> ShapeSearch<'a> {
    /// seq holds the turns made so far; pos is the end of the last edge drawn; pairs are
    /// (start state, state of the last edge) for every rim walk that fits the turns so far
    fn dfs(&mut self, turned: f64, pos: (f64, f64), pairs: &[(u32, u32)]) {
        self.nodes += 1;
        if self.nodes > self.node_cap { return; }
        let net = self.net;
        for c in self.c1..net.cls.len() {
            if self.seq.len() == 1 && c != self.second { continue; }
            let t2 = turned + net.cls[c];
            let rho = TWO_PI - t2;
            if rho < -CLOSE { break; }
            if rho.abs() < CLOSE {
                if pos.0.abs() < CLOSE && pos.1.abs() < CLOSE && self.seq.len() + 1 >= self.nmin
                    && pairs.iter().any(|&(s, cur)| net.step(cur, c as u16).any(|to| to == s)) {
                    let mut w = self.seq.clone();
                    w.push(c as u16);
                    self.found.push(w);
                }
                continue;
            }
            // every later turn is at least the first one
            if rho < self.t1 - CLOSE { continue; }
            let kmax = ((rho + CLOSE) / self.t1).floor() as usize;
            if self.seq.len() + 1 + kmax < self.nmin || self.seq.len() + 2 > self.nmax { continue; }
            let p2 = (pos.0 + t2.cos(), pos.1 + t2.sin());
            let dist = (p2.0 * p2.0 + p2.1 * p2.1).sqrt();
            if dist > kmax as f64 - 1.0 + 1e-6 { continue; }
            if rho < PI - 1e-6 && dist > 1e-6 {
                let mut ang = (-p2.1).atan2(-p2.0) - t2;
                while ang < -1e-6 { ang += TWO_PI; }
                while ang >= TWO_PI - 1e-6 { ang -= TWO_PI; }
                if ang > rho + 1e-6 { continue; }
            }
            let mut np: Vec<(u32, u32)> = Vec::new();
            for &(s, cur) in pairs {
                for to in net.step(cur, c as u16) {
                    if t2 + net.dist(to, s) <= TWO_PI + CLOSE { np.push((s, to)); }
                }
            }
            if np.is_empty() { continue; }
            np.sort_unstable();
            np.dedup();
            self.seq.push(c as u16);
            self.dfs(t2, p2, &np);
            self.seq.pop();
        }
    }
}

fn place(lib: &[Solid], g: &Graph, st: u32, pos: (f64, f64), h: f64) -> Piece {
    let s = &g.states[st as usize];
    let rz: M3 = [[h.cos(), -h.sin(), 0.0], [h.sin(), h.cos(), 0.0], [0.0, 0.0, 1.0]];
    let x = Xf { r: rz, t: [pos.0, pos.1, 0.0] }.mul(&s.xf);
    let v: Vec<V3> = lib[s.solid].v.iter().map(|p| x.ap(*p)).collect();
    let mut c = [0.0; 3];
    for p in &v { c = add(c, *p); }
    c = scale(c, 1.0 / v.len() as f64);
    Piece { solid: s.solid, state: st, r: x.r, v, c }
}

/// vertices in the window plane must lie inside the polygon; returns how many there are
fn in_plane_ok(q: &Piece, poly: &[(f64, f64)]) -> Option<usize> {
    let n = poly.len();
    let mut k = 0;
    for p in &q.v {
        if p[2] > 1e-7 { return None; }
        if p[2] > -1e-7 {
            k += 1;
            for i in 0..n {
                let (a, b) = (poly[i], poly[(i + 1) % n]);
                if (b.0 - a.0) * (p[1] - a.1) - (b.1 - a.1) * (p[0] - a.0) < -1e-7 { return None; }
            }
        }
    }
    Some(k)
}

fn same_piece(a: &Piece, b: &Piece) -> bool {
    a.solid == b.solid && dist(a.c, b.c) < 1e-6 && a.v.iter().all(|p| b.v.iter().any(|r| dist(*p, *r) < 1e-6))
}

struct Shape {
    turns: Vec<f64>,
    poly: Vec<(f64, f64)>,
    heads: Vec<f64>,
    /// per rim edge, the states that occur there in some closed rim walk
    options: Vec<Vec<u32>>,
    walks: f64,
    example: Option<(Vec<u32>, Vec<Piece>)>,
    tried: u64,
}

fn canon<T: Ord + Copy>(seq: &[T]) -> Vec<T> {
    let n = seq.len();
    (0..n).map(|k| (0..n).map(|i| seq[(i + k) % n]).collect::<Vec<T>>()).min().unwrap()
}

fn clashes(lib: &[Solid], placed: &[Piece], q: &Piece) -> bool {
    placed.iter().any(|o| !same_piece(o, q) && dist(o.c, q.c) <= lib[o.solid].radius + lib[q.solid].radius && sat(lib, o, q) < -1e-6)
}

/// the hidden piece of a corner, placed for the rim edge at pos with heading h
fn place_via(lib: &[Solid], g: &Graph, arc: u32, pos: (f64, f64), h: f64) -> Option<Piece> {
    let (si, vx) = g.arcs[arc as usize].via?;
    let rz: M3 = [[h.cos(), -h.sin(), 0.0], [h.sin(), h.cos(), 0.0], [0.0, 0.0, 1.0]];
    let x = Xf { r: rz, t: [pos.0, pos.1, 0.0] }.mul(&vx);
    let v: Vec<V3> = lib[si].v.iter().map(|p| x.ap(*p)).collect();
    let mut c = [0.0; 3];
    for p in &v { c = add(c, *p); }
    c = scale(c, 1.0 / v.len() as f64);
    Some(Piece { solid: si, state: u32::MAX, r: x.r, v, c })
}

struct Ex<'a> { lib: &'a [Solid], g: &'a Graph, net: &'a Net, seq: &'a [u16], sh: &'a Shape, lay: &'a [Vec<u32>], rim: Vec<u32>, arcs: Vec<u32>, placed: Vec<Piece>, budget: u64 }

impl<'a> Ex<'a> {
    /// depth-first search for one rim without overlapping pieces; arcs[i] is the corner after rim edge i
    fn go(&mut self, i: usize) -> bool {
        let n = self.seq.len();
        let none = Vec::new();
        if i == n {
            // the closing corner
            let cands = self.net.between.get(&(self.rim[n - 1], self.rim[0], self.seq[n - 1])).unwrap_or(&none).clone();
            for arc in cands {
                match place_via(self.lib, self.g, arc, self.sh.poly[n - 1], self.sh.heads[n - 1]) {
                    None => { self.arcs.push(arc); return true; }
                    Some(v) => if in_plane_ok(&v, &self.sh.poly).is_some() && !clashes(self.lib, &self.placed, &v) { self.arcs.push(arc); self.placed.push(v); return true; }
                }
            }
            return false;
        }
        for &st in &self.lay[i] {
            let cands: Vec<Option<u32>> = if i == 0 { vec![None] } else { self.net.between.get(&(self.rim[i - 1], st, self.seq[i - 1])).unwrap_or(&none).iter().map(|x| Some(*x)).collect() };
            if cands.is_empty() { continue; }
            if self.budget == 0 { return false; }
            self.budget -= 1;
            let q = place(self.lib, self.g, st, self.sh.poly[i], self.sh.heads[i]);
            if clashes(self.lib, &self.placed, &q) { continue; }
            for arc in cands {
                let mark = self.placed.len();
                if let Some(a) = arc {
                    if let Some(v) = place_via(self.lib, self.g, a, self.sh.poly[i - 1], self.sh.heads[i - 1]) {
                        if in_plane_ok(&v, &self.sh.poly).is_none() || clashes(self.lib, &self.placed, &v) || (!same_piece(&v, &q) && sat(self.lib, &v, &q) < -1e-6) { continue; }
                        self.placed.push(v);
                    }
                    self.arcs.push(a);
                }
                self.rim.push(st);
                self.placed.push(place(self.lib, self.g, st, self.sh.poly[i], self.sh.heads[i]));
                if self.go(i + 1) { return true; }
                self.placed.truncate(mark);
                self.rim.pop();
                if arc.is_some() { self.arcs.pop(); }
            }
        }
        false
    }
}

fn analyse(lib: &[Solid], g: &Graph, net: &Net, seq: &[u16], budget: u64) -> Shape {
    let n = seq.len();
    let turns: Vec<f64> = seq.iter().map(|c| net.cls[*c as usize]).collect();
    let (mut poly, mut heads) = (Vec::new(), Vec::new());
    let (mut h, mut pos) = (0.0f64, (0.0f64, 0.0f64));
    for i in 0..n { poly.push(pos); heads.push(h); pos = (pos.0 + h.cos(), pos.1 + h.sin()); h += turns[i]; }
    let mut okc: Vec<HashMap<u32, bool>> = vec![HashMap::new(); n];
    let mut fits = |i: usize, st: u32| -> bool { *okc[i].entry(st).or_insert_with(|| in_plane_ok(&place(lib, g, st, poly[i], heads[i]), &poly).is_some()) };
    let mut options: Vec<std::collections::BTreeSet<u32>> = vec![Default::default(); n];
    let mut walks = 0.0;
    let mut per_start: Vec<(u32, Vec<Vec<u32>>)> = Vec::new();
    for s in 0..net.ns as u32 {
        if net.cidx[s as usize] == u32::MAX || net.step(s, seq[0]).next().is_none() || !fits(0, s) { continue; }
        let mut f: Vec<Vec<u32>> = vec![vec![s]];
        for i in 0..n {
            let mut nx: Vec<u32> = f[i].iter().flat_map(|&x| net.step(x, seq[i])).collect();
            nx.sort_unstable();
            nx.dedup();
            if i + 1 < n { nx.retain(|&x| fits(i + 1, x)); } else { nx.retain(|&x| x == s); }
            if nx.is_empty() { break; }
            f.push(nx);
        }
        if f.len() != n + 1 { continue; }
        // backward: keep what can still reach s, and count walks
        let mut cnt: HashMap<u32, f64> = HashMap::new();
        cnt.insert(s, 1.0);
        let mut lay: Vec<Vec<u32>> = vec![Vec::new(); n];
        for i in (0..n).rev() {
            let mut nc: HashMap<u32, f64> = HashMap::new();
            for &x in &f[i] {
                let c: f64 = net.step(x, seq[i]).filter_map(|t| cnt.get(&t)).sum();
                if c > 0.0 { nc.insert(x, c); lay[i].push(x); options[i].insert(x); }
            }
            cnt = nc;
        }
        walks += cnt.get(&s).copied().unwrap_or(0.0);
        // the example should use the plainest solids that work
        for l in lay.iter_mut() { l.sort_by_key(|x| (lib[g.states[*x as usize].solid].faces.len(), *x)); }
        per_start.push((s, lay));
    }
    let mut sh = Shape { turns, poly, heads, options: options.into_iter().map(|o| o.into_iter().collect()).collect(), walks, example: None, tried: 0 };
    let mut left = budget;
    for (_, lay) in &per_start {
        let mut ex = Ex { lib, g, net, seq, sh: &sh, lay, rim: Vec::new(), arcs: Vec::new(), placed: Vec::new(), budget: left };
        let hit = ex.go(0);
        left = ex.budget;
        let (rim, placed) = (ex.rim, ex.placed);
        if hit { sh.example = Some((rim, placed)); break; }
        if left == 0 { break; }
    }
    sh.tried = budget - left;
    sh
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |k: &str| args.iter().position(|a| a == k).map(|i| args[i + 1].clone());
    let has = |k: &str| args.iter().any(|a| a == k);
    let lib = solids::load(&get("--solids").unwrap_or("solids.txt".into()));
    let mut active = vec![true; lib.len()];
    if let Some(o) = get("--only") {
        let names: Vec<&str> = o.split(',').collect();
        for (i, s) in lib.iter().enumerate() { active[i] = names.contains(&s.name.as_str()); }
    }
    if let Some(m) = get("--max-faces") { let m: usize = m.parse().unwrap(); for (i, s) in lib.iter().enumerate() { if s.faces.len() > m { active[i] = false; } } }
    let out = get("--out").unwrap_or("windows/out".into());
    let nmin: usize = get("--min-n").map(|x| x.parse().unwrap()).unwrap_or(13);
    let nmax: usize = get("--max-n").map(|x| x.parse().unwrap()).unwrap_or(64);
    let node_cap: u64 = get("--node-cap").map(|x| x.parse().unwrap()).unwrap_or(u64::MAX);
    let budget: u64 = get("--example-budget").map(|x| x.parse().unwrap()).unwrap_or(200_000);
    let threads: usize = get("--threads").map(|x| x.parse().unwrap()).unwrap_or(std::thread::available_parallelism().map(|x| x.get()).unwrap_or(1));
    std::fs::create_dir_all(&out).unwrap();
    let t0 = std::time::Instant::now();
    let g = build(&lib, &active, !has("--no-via"));
    let net = make_net(&g, has("--allow-coplanar"), has("--no-faces-in-plane"));
    eprintln!("states {} corners {}; usable states {} corners {} distinct turns {} smallest {:.4} deg; {:.1}s",
        g.states.len(), g.arcs.len(), net.na, net.out.iter().map(|x| x.len()).sum::<usize>(), net.cls.len(), net.cls.first().map_or(0.0, |x| x.to_degrees()), t0.elapsed().as_secs_f64());
    // one job per (first turn, second turn); the first turn is the smallest of the polygon
    let all_pairs: Vec<(u32, u32)> = (0..net.ns as u32).filter(|s| net.cidx[*s as usize] != u32::MAX).map(|s| (s, s)).collect();
    let mut jobs: Vec<(usize, usize)> = Vec::new();
    for c1 in (0..net.cls.len()).filter(|c| net.cls[*c] * nmin as f64 <= TWO_PI + CLOSE) { for c2 in c1..net.cls.len() { jobs.push((c1, c2)); } }
    let done = std::sync::atomic::AtomicUsize::new(0);
    let next = std::sync::atomic::AtomicUsize::new(0);
    let mut seqs: Vec<Vec<u16>> = Vec::new();
    let (mut nodes, mut capped) = (0u64, false);
    std::thread::scope(|sc| {
        let hs: Vec<_> = (0..threads).map(|_| sc.spawn(|| {
            let (mut mine, mut nd, mut cp) = (Vec::new(), 0u64, false);
            loop {
                let j = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if j >= jobs.len() { break; }
                let (c1, c2) = jobs[j];
                let mut se = ShapeSearch { net: &net, nmin, nmax, t1: net.cls[c1], c1, seq: Vec::new(), found: Vec::new(), nodes: 0, node_cap, second: c2 };
                // the first turn: handled by the general step with only class c1 allowed
                let t2 = net.cls[c1];
                let mut np: Vec<(u32, u32)> = Vec::new();
                for &(s, cur) in &all_pairs { for to in net.step(cur, c1 as u16) { if t2 + net.dist(to, s) <= TWO_PI + CLOSE { np.push((s, to)); } } }
                np.sort_unstable();
                np.dedup();
                if !np.is_empty() {
                    se.seq.push(c1 as u16);
                    se.dfs(t2, (1.0 + t2.cos(), t2.sin()), &np);
                }
                let dn = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                if se.nodes > 200_000 || dn % 20_000 == 0 { eprintln!("  turns {:.4} {:.4}: {} closed turn sequences, {} nodes; {} of {} jobs, {:.0}s", net.cls[c1].to_degrees(), net.cls[c2].to_degrees(), se.found.len(), se.nodes, dn, jobs.len(), t0.elapsed().as_secs_f64()); }
                nd += se.nodes;
                cp |= se.nodes > node_cap;
                mine.append(&mut se.found);
            }
            (mine, nd, cp)
        })).collect();
        for h in hs { let (mut m, nd, cp) = h.join().unwrap(); seqs.append(&mut m); nodes += nd; capped |= cp; }
    });
    let mut seen = std::collections::HashSet::new();
    seqs.retain(|w| seen.insert(canon(w)));
    seqs.sort_by_key(|w| (w.len(), w.clone()));
    eprintln!("turn sequences {} (one per rotation) nodes {} in {:.1}s{}", seqs.len(), nodes, t0.elapsed().as_secs_f64(), if capped { " (NODE CAP HIT, incomplete)" } else { "" });
    // analyse each shape
    let nexts = std::sync::atomic::AtomicUsize::new(0);
    let mut shapes: Vec<(usize, Shape)> = Vec::new();
    std::thread::scope(|sc| {
        let hs: Vec<_> = (0..threads).map(|_| sc.spawn(|| {
            let mut mine = Vec::new();
            loop {
                let j = nexts.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if j >= seqs.len() { break; }
                mine.push((j, analyse(&lib, &g, &net, &seqs[j], budget)));
            }
            mine
        })).collect();
        for h in hs { shapes.append(&mut h.join().unwrap()); }
    });
    shapes.sort_by_key(|x| x.0);
    let mut outline: HashMap<Vec<u16>, usize> = HashMap::new();
    let mut f = std::io::BufWriter::new(std::fs::File::create(format!("{}/windows.jsonl", out)).unwrap());
    let mut byn: std::collections::BTreeMap<usize, [usize; 4]> = Default::default();
    let mut id = 0;
    for (j, sh) in &shapes {
        let e = byn.entry(seqs[*j].len()).or_default();
        if sh.walks == 0.0 { e[2] += 1; continue; }
        if sh.example.is_none() { e[3] += 1; }
        let seq = &seqs[*j];
        let n = seq.len();
        let mut rv = seq.clone();
        rv.reverse();
        let no = outline.len();
        let ol = *outline.entry(canon(seq).min(canon(&rv))).or_insert(no);
        if ol == no { e[1] += 1; }
        e[0] += 1;
        let rot = (1..=n).filter(|k| (0..n).all(|i| seq[i] == seq[(i + k) % n])).count();
        let mirror = (0..n).any(|k| (0..n).all(|i| seq[i] == rv[(i + k) % n]));
        let jl = |v: Vec<String>| v.join(",");
        let ex = match &sh.example {
            None => "null".to_string(),
            Some((rim, placed)) => {
                let mut pcs: Vec<&Piece> = Vec::new();
                for q in placed { if !pcs.iter().any(|o| same_piece(o, q)) { pcs.push(q); } }
                let ep: Vec<String> = (0..n).map(|i| { let q = place(&lib, &g, rim[i], sh.poly[i], sh.heads[i]); pcs.iter().position(|o| same_piece(o, &q)).unwrap().to_string() }).collect();
                let fip = pcs.iter().filter(|q| in_plane_ok(q, &sh.poly).unwrap_or(0) > 2).count();
                let hidden = pcs.iter().filter(|q| q.state == u32::MAX).count();
                let ps = pcs.iter().map(|p| {
                    let t = sub(p.v[0], mv(&p.r, lib[p.solid].v[0]));
                    let x: Vec<String> = p.r.iter().flatten().chain(t.iter()).map(|v| format!("{:.12}", v)).collect();
                    format!("{{\"solid\":\"{}\",\"hidden\":{},\"xf\":[{}]}}", lib[p.solid].name, p.state == u32::MAX, x.join(","))
                }).collect();
                format!("{{\"states\":[{}],\"edge_piece\":[{}],\"faces_in_plane\":{},\"hidden\":{},\"pieces\":[{}]}}", jl(rim.iter().map(|x| x.to_string()).collect()), jl(ep), fip, hidden, jl(ps))
            }
        };
        let rmin = {
            // inradius about the centroid: a rough size for sorting
            let c = sh.poly.iter().fold((0.0, 0.0), |a, p| (a.0 + p.0 / n as f64, a.1 + p.1 / n as f64));
            (0..n).map(|i| { let (a, h) = (sh.poly[i], sh.heads[i]); -((c.0 - a.0) * h.sin() - (c.1 - a.1) * h.cos()) }).fold(f64::MAX, f64::min).abs()
        };
        writeln!(f, "{{\"id\":{},\"n\":{},\"outline\":{},\"rot\":{},\"mirror\":{},\"inradius\":{:.6},\"rims\":{:e},\"verified\":{},\"tried\":{},\"turns\":[{}],\"poly\":[{}],\"options\":[{}],\"example\":{}}}",
            id, n, ol, rot, mirror, rmin, sh.walks, sh.example.is_some(), sh.tried,
            jl(sh.turns.iter().map(|t| format!("{:.4}", t.to_degrees())).collect()),
            jl(sh.poly.iter().map(|p| format!("[{:.9},{:.9}]", p.0, p.1)).collect()),
            jl(sh.options.iter().map(|o| format!("[{}]", jl(o.iter().map(|x| x.to_string()).collect()))).collect()), ex).unwrap();
        id += 1;
    }
    // tables of the states and corners that occur in some window
    let full = has("--full-tables");
    let mut used: Vec<bool> = (0..g.states.len()).map(|i| full && net.cidx[i] != u32::MAX).collect();
    for (_, sh) in &shapes { if sh.walks > 0.0 { for o in &sh.options { for x in o { used[*x as usize] = true; } } } }
    drop(f);
    {
        let mut f = std::io::BufWriter::new(std::fs::File::create(format!("{}/states.tsv", out)).unwrap());
        writeln!(f, "id\tsolid\ta\tb\ttilt_deg\tspan_deg\txf").unwrap();
        for (i, s) in g.states.iter().enumerate() {
            if !used[i] { continue; }
            let x: Vec<String> = s.xf.r.iter().flatten().chain(s.xf.t.iter()).map(|v| format!("{:.12}", v)).collect();
            writeln!(f, "{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}", i, lib[s.solid].name, s.a, s.b, s.tilt.to_degrees(), s.span.to_degrees(), x.join(" ")).unwrap();
        }
        let mut f = std::io::BufWriter::new(std::fs::File::create(format!("{}/corners.tsv", out)).unwrap());
        writeln!(f, "from\tto\tturn_deg\tglue_face\taplanar\tvia\tvia_xf").unwrap();
        for a in &g.arcs {
            if !used[a.from as usize] || !used[a.to as usize] { continue; }
            let vx = a.via.map_or("-".to_string(), |v| v.1.r.iter().flatten().chain(v.1.t.iter()).map(|x| format!("{:.12}", x)).collect::<Vec<_>>().join(" "));
            writeln!(f, "{}\t{}\t{:.6}\t{}\t{}\t{}\t{}", a.from, a.to, a.turn.to_degrees(), a.face, a.aplanar as u8, a.via.map_or("-", |v| lib[v.0].name.as_str()), vx).unwrap();
        }
    }
    eprintln!("windows written {} in {:.1}s", id, t0.elapsed().as_secs_f64());
    for (n, c) in &byn { eprintln!("  n = {:2}: {} polygons ({} up to mirror image), {} more dropped because a face sticks out of the polygon, {} without a collision-free rim found", n, c[0], c[1], c[2], c[3]); }
}
