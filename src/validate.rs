use crate::geom::*;
use crate::solids::Solid;
use std::collections::HashMap as StdMap;
use std::hash::{BuildHasherDefault, Hasher};
#[derive(Default)]
pub struct Fx(u64);
impl Hasher for Fx {
    fn finish(&self) -> u64 { self.0 ^ (self.0 >> 32) }
    fn write(&mut self, b: &[u8]) { for &x in b { self.write_u64(x as u64); } }
    fn write_u32(&mut self, x: u32) { self.write_u64(x as u64); }
    fn write_i64(&mut self, x: i64) { self.write_u64(x as u64); }
    fn write_u64(&mut self, x: u64) { self.0 = (self.0.rotate_left(5) ^ x).wrapping_mul(0x517CC1B727220A95); }
}
type HashMap<K, V> = StdMap<K, V, BuildHasherDefault<Fx>>;
use std::sync::atomic::{AtomicU64, Ordering};
pub static T: [AtomicU64; 6] = [AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)];
fn tick(i: usize, t: &mut std::time::Instant) { T[i].fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed); *t = std::time::Instant::now(); }

pub struct Placed {
    pub solid: usize,
    pub xf: Xf,
    pub v: Vec<V3>,
    pub fnorm: Vec<V3>,
    pub foff: Vec<f64>,
    pub ed: Vec<V3>,
    pub c: V3,
    pub r: f64,
    pub inr: f64,
}

pub fn place(lib: &[Solid], solid: usize, xf: Xf) -> Placed {
    let s = &lib[solid];
    let v: Vec<V3> = s.v.iter().map(|p| xf.ap(*p)).collect();
    let fnorm: Vec<V3> = s.fnorm.iter().map(|n| mv(&xf.r, *n)).collect();
    let foff: Vec<f64> = (0..s.faces.len()).map(|f| dot(fnorm[f], v[s.faces[f][0]])).collect();
    let ed: Vec<V3> = s.edirs.iter().map(|d| mv(&xf.r, *d)).collect();
    Placed { solid, xf, v, fnorm, foff, ed, c: xf.t, r: s.radius, inr: s.inr }
}

fn gap_on(a: &Placed, b: &Placed, ax: V3) -> f64 {
    let (mut amin, mut amax, mut bmin, mut bmax) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for p in &a.v { let d = dot(*p, ax); if d < amin { amin = d; } if d > amax { amax = d; } }
    for p in &b.v { let d = dot(*p, ax); if d < bmin { bmin = d; } if d > bmax { bmax = d; } }
    (bmin - amax).max(amin - bmax)
}

/// Separation test. Returns (gap, axis): gap >= -EPS means the axis separates (touching if |gap| <= EPS);
/// gap < -EPS means the interiors overlap.
pub fn sat(a: &Placed, b: &Placed) -> (f64, V3) {
    let z = [0.0; 3];
    let cc = sub(b.c, a.c);
    let d = norm(cc);
    if d < a.inr + b.inr - 1e-4 { return (-1.0, z); }
    // a centre or vertex strictly inside the other piece certifies overlap without the edge-pair axes
    for (p, q) in [(a, b), (b, a)] {
        if (0..q.fnorm.len()).all(|f| dot(q.fnorm[f], p.c) - q.foff[f] < -1e-4) { return (-1.0, z); }
        for v in &p.v {
            if (0..q.fnorm.len()).all(|f| dot(q.fnorm[f], *v) - q.foff[f] < -1e-4) { return (-1.0, z); }
        }
    }
    let mut best = f64::MIN;
    if d > 1e-9 {
        let ax = scale(cc, 1.0 / d);
        let g = gap_on(a, b, ax);
        if g >= -EPS { return (g, ax); }
    }
    for n in a.fnorm.iter().chain(b.fnorm.iter()) {
        let g = gap_on(a, b, *n);
        if g >= -EPS { return (g, *n); }
        if g > best { best = g; }
    }
    for x in &a.ed {
        for y in &b.ed {
            let c = cross(*x, *y);
            let l = norm(c);
            if l < 1e-6 { continue; }
            let ax = scale(c, 1.0 / l);
            let g = gap_on(a, b, ax);
            if g >= -EPS { return (g, ax); }
            if g > best { best = g; }
        }
    }
    (best, z)
}

/// Touching pieces separated by plane normal n: true if the contact is anything other than whole
/// shared faces/edges/vertices. Only vertices and edges lying in the contact plane can take part.
fn improper(lib: &[Solid], a: &Placed, b: &Placed, n: V3) -> bool {
    let pa: Vec<f64> = a.v.iter().map(|v| dot(*v, n)).collect();
    let pb: Vec<f64> = b.v.iter().map(|v| dot(*v, n)).collect();
    let (amin, amax) = pa.iter().fold((f64::MAX, f64::MIN), |m, x| (m.0.min(*x), m.1.max(*x)));
    let (bmin, bmax) = pb.iter().fold((f64::MAX, f64::MIN), |m, x| (m.0.min(*x), m.1.max(*x)));
    let (ona, onb): (Vec<bool>, Vec<bool>) = if bmin - amax > amin - bmax {
        (pa.iter().map(|x| *x > amax - 1e-5).collect(), pb.iter().map(|x| *x < bmin + 1e-5).collect())
    } else {
        (pa.iter().map(|x| *x < amin + 1e-5).collect(), pb.iter().map(|x| *x > bmax - 1e-5).collect())
    };
    for (p, q, on) in [(a, b, &ona), (b, a, &onb)] {
        for (vi, v) in p.v.iter().enumerate() {
            if !on[vi] { continue; }
            let mut m = f64::MIN;
            for f in 0..q.fnorm.len() { let d = dot(q.fnorm[f], *v) - q.foff[f]; if d > m { m = d; } }
            if m.abs() <= 1e-6 && !q.v.iter().any(|w| dist(*w, *v) < 1e-5) { return true; }
        }
    }
    let eb: Vec<(usize, usize)> = lib[b.solid].edges.iter().cloned().filter(|&(k, l)| onb[k] && onb[l]).collect();
    if eb.is_empty() { return false; }
    for &(i, j) in &lib[a.solid].edges {
        if !(ona[i] && ona[j]) { continue; }
        let (p1, d1) = (a.v[i], sub(a.v[j], a.v[i]));
        for &(k, l) in &eb {
            let (p2, d2) = (b.v[k], sub(b.v[l], b.v[k]));
            let n = cross(d1, d2);
            let n2 = dot(n, n);
            if n2 < 1e-10 { continue; }
            let w = sub(p2, p1);
            if (dot(w, n) / n2.sqrt()).abs() > 1e-6 { continue; }
            let s = dot(cross(w, d2), n) / n2;
            let t = dot(cross(w, d1), n) / n2;
            if s > 1e-5 && s < 1.0 - 1e-5 && t > 1e-5 && t < 1.0 - 1e-5 { return true; }
        }
    }
    false
}

/// Supporting planes of the point set that contain the segment u-v (none if the segment is not on
/// the hull). Planes are (outward unit normal, offset). All points must project into a half-plane
/// through the segment; the extreme directions give the planes.
pub fn support_planes(pts: &[V3], u: V3, v: V3, centre: V3) -> Vec<(V3, f64)> {
    let d = unit(sub(v, u));
    let t = if d[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    let e1 = unit(cross(d, t));
    let e2 = cross(d, e1);
    let mut rf: Option<(f64, f64)> = None;
    let (mut hi, mut lo): (Option<(f64, f64)>, Option<(f64, f64)>) = (None, None);
    let (mut xhi, mut xlo, mut xrf) = ([0.0; 3], [0.0; 3], [0.0; 3]);
    let mut opp = false;
    let cr = |a: (f64, f64), b: (f64, f64)| a.0 * b.1 - a.1 * b.0;
    let n2 = |a: (f64, f64)| a.0 * a.0 + a.1 * a.1;
    for p in pts {
        let x = sub(*p, u);
        let q = (dot(x, e1), dot(x, e2));
        if n2(q) < 1e-12 { continue; }
        let r = match rf { None => { rf = Some(q); xrf = x; continue; } Some(r) => r };
        let s = cr(r, q);
        if s * s <= 1e-13 * n2(r) * n2(q) {
            if r.0 * q.0 + r.1 * q.1 < 0.0 { opp = true; if hi.is_some() && lo.is_some() { return Vec::new(); } }
            continue;
        }
        if s > 0.0 {
            if hi.map_or(true, |h| cr(h, q) > 0.0) { hi = Some(q); xhi = x; }
        } else if lo.map_or(true, |l| cr(l, q) < 0.0) { lo = Some(q); xlo = x; }
        if let (Some(h), Some(l)) = (hi, lo) {
            if opp { return Vec::new(); }
            let c = cr(l, h);
            if c < 0.0 && c * c > 1e-13 * n2(l) * n2(h) { return Vec::new(); }
        }
    }
    if rf.is_none() { return Vec::new(); }
    let cands: Vec<V3> = if opp { vec![xrf] } else { vec![if lo.is_some() { xlo } else { xrf }, if hi.is_some() { xhi } else { xrf }] };
    let mut out: Vec<(V3, f64)> = Vec::new();
    for q in cands {
        let mut n = unit(cross(d, q));
        if dot(n, sub(centre, u)) > 0.0 { n = scale(n, -1.0); }
        let c = dot(n, u);
        if !out.iter().any(|(m, k)| dist(*m, n) < 1e-6 && (k - c).abs() < 1e-6) { out.push((n, c)); }
    }
    out
}

pub struct Welder { map: HashMap<(i64, i64, i64), u32>, pub pts: Vec<V3> }
impl Welder {
    pub fn new() -> Welder { Welder { map: HashMap::default(), pts: Vec::new() } }
    pub fn id(&mut self, p: V3) -> u32 {
        let mut alt = [[0i64; 2]; 3];
        let mut na = [1usize; 3];
        for i in 0..3 {
            let x = p[i] * 1e4;
            let a = x.floor();
            alt[i][0] = a as i64;
            let fr = x - a;
            if fr < 0.02 { alt[i][1] = a as i64 - 1; na[i] = 2; } else if fr > 0.98 { alt[i][1] = a as i64 + 1; na[i] = 2; }
        }
        for i in 0..na[0] { for j in 0..na[1] { for k in 0..na[2] {
            if let Some(&id) = self.map.get(&(alt[0][i], alt[1][j], alt[2][k])) { return id; }
        } } }
        let id = self.pts.len() as u32;
        self.map.insert((alt[0][0], alt[1][0], alt[2][0]), id);
        self.pts.push(p);
        id
    }
}

pub struct Res {
    pub gorder: usize,
    pub npieces: usize,
    pub counts: Vec<(usize, usize)>,
    pub genus: i64,
    pub nv: usize,
    pub ne: usize,
    pub nf: usize,
    pub a_viol: usize,
    pub q_ok: bool,
    pub q_missing: usize,
    pub q_edges: usize,
    pub hull_faces: Vec<usize>,
    pub hull_regular: bool,
    pub rin: f64,
    pub rcirc: f64,
    pub pts: Vec<V3>,
    pub faces: Vec<Vec<u32>>,
    pub fsolid: Vec<usize>,
    pub reason: &'static str,
}

/// Build the orbit of reps under gens, check it is a valid embedded closed surface, and measure it.
/// Err(reason) on rejection.
pub fn validate(lib: &[Solid], reps: &[(usize, Xf)], grp: &[Xf], centre: V3, prescreen: bool) -> Result<Res, &'static str> {
    let mut tm = std::time::Instant::now();
    // the centre strictly inside an off-centre piece means that piece meets its own images
    for &(s, x) in reps {
        let o = mv(&mt(&x.r), sub(centre, x.t));
        let sl = &lib[s];
        if norm(o) > 1e-6 && (0..sl.faces.len()).all(|f| dot(sl.fnorm[f], sub(o, sl.v[sl.faces[f][0]])) < -1e-6) { return Err("overlap-centre"); }
    }
    tick(0, &mut tm);
    // collisions: representatives against each other and against every image, built lazily
    let rp: Vec<Placed> = reps.iter().map(|&(s, x)| place(lib, s, x)).collect();
    let mut touching: Vec<(usize, Placed, V3)> = Vec::new();
    let mut rtouch: Vec<(usize, usize, V3)> = Vec::new();
    for i in 0..rp.len() {
        for j in 0..i {
            let d = dist(rp[i].c, rp[j].c);
            if d < 1e-5 { return Err("overlap-chain"); }
            if d <= rp[i].r + rp[j].r + 1e-6 {
                let (g, ax) = sat(&rp[i], &rp[j]);
                if g < -EPS { return Err("overlap-chain"); }
                if g <= EPS { rtouch.push((i, j, ax)); }
            }
        }
    }
    // cheapest certificate first: inscribed spheres of a representative and an image overlap
    for g in &grp[1..] {
        for &(s, x) in reps {
            let c = g.ap(x.t);
            for a in &rp {
                let d = dist(a.c, c);
                if d > 1e-5 && d < a.inr + lib[s].inr - 1e-4 { return Err("overlap-sphere"); }
            }
        }
    }
    for g in &grp[1..] {
        for &(s, x) in reps {
            let xf = g.mul(&x);
            let mut img: Option<Placed> = None;
            for (ai, a) in rp.iter().enumerate() {
                let d = dist(a.c, xf.t);
                if d > a.r + lib[s].radius + 1e-6 { continue; }
                if d < 1e-5 {
                    if a.solid != s { return Err("overlap"); }
                    let rel = mmul(&mt(&a.xf.r), &xf.r);
                    if !lib[s].syms.iter().any(|m| mclose(m, &rel, 1e-6)) { return Err("overlap"); }
                } else {
                    if img.is_none() { img = Some(place(lib, s, xf)); }
                    let (gp, ax) = sat(a, img.as_ref().unwrap());
                    if gp < -EPS { return Err("overlap"); }
                    if gp <= EPS { touching.push((ai, place(lib, s, xf), ax)); }
                }
            }
        }
    }
    for (i, j, ax) in &rtouch { if improper(lib, &rp[*i], &rp[*j], *ax) { return Err("contact"); } }
    for (i, im, ax) in &touching { if improper(lib, &rp[*i], im, *ax) { return Err("contact"); } }
    tick(1, &mut tm);
    // orbit as bare points (one flat buffer); pieces are only built if the pre-screens pass
    let nrep = rp.len();
    let mut cs: Vec<(usize, Xf)> = reps.to_vec();
    for g in &grp[1..] {
        for &(s, x) in reps {
            let xf = g.mul(&x);
            if cs.iter().any(|p| dist(p.1.t, xf.t) < 1e-5) { continue; }
            cs.push((s, xf));
        }
    }
    let mut all: Vec<V3> = Vec::with_capacity(cs.iter().map(|c| lib[c.0].v.len()).sum());
    for &(s, xf) in &cs { for v in &lib[s].v { all.push(xf.ap(*v)); } }
    // (Q) pre-screen at the outermost vertex, which is certainly a hull vertex: its unit-distance
    // neighbours must span a cone containing every other vertex, or some hull edge there is not unit.
    if prescreen {
        let mut vs = centre;
        let mut far = -1.0;
        for v in &all { let d = dist(*v, centre); if d > far { far = d; vs = *v; } }
        let mut dirs: Vec<V3> = Vec::new();
        for v in &all {
            let d = sub(*v, vs);
            if (dot(d, d) - 1.0).abs() < 2e-6 && !dirs.iter().any(|x| dist(*x, d) < 1e-6) { dirs.push(d); }
        }
        let mut facets: Vec<V3> = Vec::new();
        for i in 0..dirs.len() {
            for j in 0..i {
                let n = cross(dirs[i], dirs[j]);
                if norm(n) < 1e-7 { continue; }
                let n = unit(n);
                let (mut pos, mut neg) = (false, false);
                for d in &dirs { let x = dot(*d, n); if x > 1e-7 { pos = true; } if x < -1e-7 { neg = true; } }
                if pos && neg { continue; }
                let n = if pos { scale(n, -1.0) } else { n };
                if !facets.iter().any(|f| dist(*f, n) < 1e-6) { facets.push(n); }
            }
        }
        if facets.len() < 3 { return Err("q-prescreen"); }
        for v in &all {
            let d = sub(*v, vs);
            if facets.iter().any(|n| dot(d, *n) > 1e-6) { return Err("q-prescreen"); }
        }
    }
    // (Q) necessary condition on bare points, no mesh: every hull face found through a
    // representative edge must have unit-length sides. Stops at the first failure.
    if prescreen {
        let mut done: Vec<(V3, f64)> = Vec::new();
        for r in &rp {
            for &(ia, ib) in &lib[r.solid].edges {
                for (n, c) in support_planes(&all, r.v[ia], r.v[ib], centre) {
                    if done.iter().any(|(m, k)| dist(*m, n) < 1e-6 && (k - c).abs() < 1e-6) { continue; }
                    done.push((n, c));
                    if all.iter().any(|p| dot(n, *p) - c > 1e-6) { continue; }
                    let t = if n[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
                    let e1 = unit(cross(n, t));
                    let e2 = cross(n, e1);
                    let mut on: Vec<(f64, f64)> = all.iter().filter(|p| (dot(n, **p) - c).abs() < 1e-6).map(|p| (dot(*p, e1), dot(*p, e2))).collect();
                    on.sort_by(|a, b| ((a.0 * 1e6).round(), (a.1 * 1e6).round()).partial_cmp(&((b.0 * 1e6).round(), (b.1 * 1e6).round())).unwrap());
                    let cr = |o: &(f64, f64), a: &(f64, f64), b: &(f64, f64)| (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0);
                    let mut hull: Vec<(f64, f64)> = Vec::new();
                    for pass in 0..2 {
                        let start = hull.len();
                        for k in 0..on.len() {
                            let p = if pass == 0 { on[k] } else { on[on.len() - 1 - k] };
                            while hull.len() >= start + 2 && cr(&hull[hull.len() - 2], &hull[hull.len() - 1], &p) <= 1e-9 { hull.pop(); }
                            hull.push(p);
                        }
                        hull.pop();
                    }
                    for i in 0..hull.len() {
                        let (a, b) = (hull[i], hull[(i + 1) % hull.len()]);
                        let l = ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt();
                        // coincident points from neighbouring pieces can survive the hull as zero-length sides
                        if l > 1e-5 && (l - 1.0).abs() > 1e-5 { return Err("q-proxy"); }
                    }
                }
            }
        }
    }
    // full pieces for the mesh
    let mut pcs: Vec<Placed> = rp;
    for &(s, xf) in &cs[nrep..] { pcs.push(place(lib, s, xf)); }
    // boundary mesh
    let mut w = Welder::new();
    let mut fmap: HashMap<(u32, u64, u64), (u8, Vec<u32>, usize)> = HashMap::default();
    let mut rep_edges: Vec<(u32, u32)> = Vec::new();
    for (pi, p) in pcs.iter().enumerate() {
        let ids: Vec<u32> = p.v.iter().map(|v| w.id(*v)).collect();
        for f in &lib[p.solid].faces {
            let fi: Vec<u32> = f.iter().map(|&i| ids[i]).collect();
            // order-independent fingerprint of the vertex set
            let mut key = (u32::MAX, 0u64, 0u64);
            for &x in &fi {
                key.0 = key.0.min(x);
                let h = (x as u64 + 1).wrapping_mul(0x9E3779B97F4A7C15);
                key.1 = key.1.wrapping_add(h ^ (h >> 29));
                key.2 ^= h.rotate_left(23).wrapping_mul(0xD6E8FEB86659FD93);
            }
            match fmap.get_mut(&key) {
                Some(e) => e.0 += 1,
                None => { fmap.insert(key, (1, fi, p.solid)); }
            }
        }
        if pi < nrep {
            for &(a, b) in &lib[p.solid].edges { rep_edges.push((ids[a], ids[b])); }
        }
    }
    let mut faces: Vec<Vec<u32>> = Vec::new();
    let mut fsolid = Vec::new();
    let mut keys: Vec<&(u32, u64, u64)> = fmap.keys().collect();
    keys.sort();
    for k in keys {
        let (c, f, s) = &fmap[k];
        if *c > 2 { return Err("face3"); }
        if *c == 1 { faces.push(f.clone()); fsolid.push(*s); }
    }
    if faces.is_empty() { return Err("empty"); }
    let pts = w.pts;
    let mut de: HashMap<(u32, u32), u32> = HashMap::default();
    for (fi, f) in faces.iter().enumerate() {
        for i in 0..f.len() {
            if de.insert((f[i], f[(i + 1) % f.len()]), fi as u32).is_some() { return Err("nonmanifold"); }
        }
    }
    for (&(a, b), _) in &de { if !de.contains_key(&(b, a)) { return Err("nonmanifold"); } }
    // vertex fans
    let mut corners = vec![0u32; pts.len()];
    let mut first: Vec<Option<(u32, u32)>> = vec![None; pts.len()];
    for f in &faces {
        for i in 0..f.len() {
            let v = f[i] as usize;
            corners[v] += 1;
            if first[v].is_none() { first[v] = Some((f[(i + 1) % f.len()], f[(i + f.len() - 1) % f.len()])); }
        }
    }
    for v in 0..pts.len() {
        let (n0, mut prev) = match first[v] { Some(x) => x, None => return Err("nonmanifold") };
        let mut cnt = 1;
        // walk: next corner at v is in the face holding directed edge (v -> prev)
        while prev != n0 {
            let f = &faces[de[&(v as u32, prev)] as usize];
            let i = f.iter().position(|&x| x == v as u32).unwrap();
            prev = f[(i + f.len() - 1) % f.len()];
            cnt += 1;
            if cnt > corners[v] { break; }
        }
        if cnt != corners[v] { return Err("pinch"); }
    }
    // connected
    let mut par: Vec<usize> = (0..pts.len()).collect();
    fn find(p: &mut Vec<usize>, mut x: usize) -> usize { while p[x] != x { p[x] = p[p[x]]; x = p[x]; } x }
    for (&(a, b), _) in &de { let (x, y) = (find(&mut par, a as usize), find(&mut par, b as usize)); par[x] = y; }
    let root = find(&mut par, 0);
    for v in 0..pts.len() { if find(&mut par, v) != root { return Err("disconnected"); } }
    let (nv, ne, nf) = (pts.len(), de.len() / 2, faces.len());
    let chi = nv as i64 - ne as i64 + nf as i64;
    if chi % 2 != 0 { return Err("euler"); }
    let genus = (2 - chi) / 2;
    if genus < 1 { return Err("genus0"); }
    tick(2, &mut tm);
    // (A)
    let fnorm: Vec<V3> = faces.iter().map(|f| {
        let (a, b, c) = (pts[f[0] as usize], pts[f[1] as usize], pts[f[2] as usize]);
        unit(cross(sub(b, a), sub(c, b)))
    }).collect();
    let mut a_viol = 0;
    for (&(a, b), &f1) in &de {
        if a < b {
            let f2 = de[&(b, a)];
            if dot(fnorm[f1 as usize], fnorm[f2 as usize]) > 1.0 - 1e-9 { a_viol += 1; }
        }
    }
    // (Q): supporting planes through representative edges
    let mut planes: Vec<(V3, f64)> = Vec::new();
    for &(ia, ib) in &rep_edges {
        if !de.contains_key(&(ia, ib)) { continue; }
        let (u, v) = (pts[ia as usize], pts[ib as usize]);
        for (n, c) in support_planes(&pts, u, v, centre) {
            if !planes.iter().any(|(m, k)| dist(*m, n) < 1e-6 && (k - c).abs() < 1e-6) { planes.push((n, c)); }
        }
    }
    let rcirc = pts.iter().map(|p| dist(*p, centre)).fold(0.0, f64::max);
    let mut q_missing = 0;
    let mut q_edges = 0;
    let mut hull_regular = true;
    let mut hull_faces = Vec::new();
    let mut rin = f64::MAX;
    for &(n, c) in &planes {
        // safety: must be supporting
        if pts.iter().any(|p| dot(n, *p) - c > 1e-6) { continue; }
        rin = rin.min(c - dot(n, centre));
        let t = if n[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let e1 = unit(cross(n, t));
        let e2 = cross(n, e1);
        let mut on: Vec<(f64, f64, u32)> = Vec::new();
        for (i, p) in pts.iter().enumerate() {
            if (dot(n, *p) - c).abs() < 1e-6 { on.push((dot(*p, e1), dot(*p, e2), i as u32)); }
        }
        on.sort_by(|a, b| ((a.0 * 1e6).round(), (a.1 * 1e6).round()).partial_cmp(&((b.0 * 1e6).round(), (b.1 * 1e6).round())).unwrap());
        let cr = |o: &(f64, f64, u32), a: &(f64, f64, u32), b: &(f64, f64, u32)| (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0);
        let mut hull: Vec<(f64, f64, u32)> = Vec::new();
        for pass in 0..2 {
            let start = hull.len();
            let it: Box<dyn Iterator<Item = &(f64, f64, u32)>> = if pass == 0 { Box::new(on.iter()) } else { Box::new(on.iter().rev()) };
            for p in it {
                while hull.len() >= start + 2 && cr(&hull[hull.len() - 2], &hull[hull.len() - 1], p) <= 1e-9 { hull.pop(); }
                hull.push(*p);
            }
            hull.pop();
        }
        if hull.len() < 3 { continue; }
        hull_faces.push(hull.len());
        let (cx, cy) = (hull.iter().map(|h| h.0).sum::<f64>() / hull.len() as f64, hull.iter().map(|h| h.1).sum::<f64>() / hull.len() as f64);
        let r0 = ((hull[0].0 - cx).powi(2) + (hull[0].1 - cy).powi(2)).sqrt();
        for i in 0..hull.len() {
            let (a, b) = (hull[i], hull[(i + 1) % hull.len()]);
            q_edges += 1;
            if !de.contains_key(&(a.2, b.2)) { q_missing += 1; }
            let r = ((a.0 - cx).powi(2) + (a.1 - cy).powi(2)).sqrt();
            if (r - r0).abs() > 1e-6 { hull_regular = false; }
        }
    }
    tick(3, &mut tm);
    hull_faces.sort();
    let q_ok = !hull_faces.is_empty() && q_missing == 0;
    let mut counts: Vec<(usize, usize)> = Vec::new();
    for p in &pcs {
        match counts.iter_mut().find(|c| c.0 == p.solid) { Some(c) => c.1 += 1, None => counts.push((p.solid, 1)) }
    }
    counts.sort();
    Ok(Res { gorder: grp.len(), npieces: pcs.len(), counts, genus, nv, ne, nf, a_viol, q_ok, q_missing, q_edges, hull_faces, hull_regular, rin, rcirc, pts, faces, fsolid, reason: "" })
}
