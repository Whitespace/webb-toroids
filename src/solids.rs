use crate::geom::*;

pub struct Solid {
    pub name: String,
    pub v: Vec<V3>,
    pub faces: Vec<Vec<usize>>,
    pub fnorm: Vec<V3>,
    pub edges: Vec<(usize, usize)>,
    pub edirs: Vec<V3>,
    pub radius: f64,
    pub inr: f64,
    pub syms: Vec<M3>,
    pub fperm: Vec<Vec<usize>>,
    pub axes: Vec<(V3, usize)>,
    /// frames[f][k]: canonical polygon frame -> face f with vertex k on +x, +z outward.
    pub frames: Vec<Vec<Xf>>,
    /// flags (face, k) up to rotational symmetry
    pub flags: Vec<(usize, usize)>,
}

fn find_vertex(v: &[V3], p: V3) -> Option<usize> { v.iter().position(|q| dist(*q, p) < 1e-6) }

impl Solid {
    fn new(name: String, v: Vec<V3>, faces: Vec<Vec<usize>>) -> Solid {
        let mut fnorm = Vec::new();
        let mut frames = Vec::new();
        for f in &faces {
            let n = f.len();
            let mut c = [0.0; 3];
            for &i in f { c = add(c, v[i]); }
            c = scale(c, 1.0 / n as f64);
            let mut nn = [0.0; 3];
            for i in 0..n { nn = add(nn, cross(sub(v[f[i]], c), sub(v[f[(i + 1) % n]], c))); }
            let e3 = unit(nn);
            assert!(dot(e3, c) > 0.0);
            fnorm.push(e3);
            let mut fr = Vec::new();
            for k in 0..n {
                let e1 = unit(sub(v[f[k]], c));
                let e2 = cross(e3, e1);
                fr.push(Xf { r: from_cols(e1, e2, e3), t: c });
            }
            frames.push(fr);
        }
        let mut edges = Vec::new();
        for f in &faces {
            for i in 0..f.len() {
                let (a, b) = (f[i], f[(i + 1) % f.len()]);
                if a < b { edges.push((a, b)); }
            }
        }
        let mut edirs: Vec<V3> = Vec::new();
        for &(a, b) in &edges {
            let d = unit(sub(v[b], v[a]));
            if !edirs.iter().any(|e| dot(*e, d).abs() > 1.0 - 1e-9) { edirs.push(d); }
        }
        let inr = (0..faces.len()).map(|f| dot(fnorm[f], v[faces[f][0]])).fold(f64::MAX, f64::min);
        let radius = v.iter().map(|p| norm(*p)).fold(0.0, f64::max);
        // rotational symmetries
        let mut syms: Vec<M3> = Vec::new();
        let mut vperms: Vec<Vec<usize>> = Vec::new();
        let n0 = faces[0].len();
        let base = mt(&frames[0][0].r);
        for f in 0..faces.len() {
            if faces[f].len() != n0 { continue; }
            for k in 0..n0 {
                let r = mmul(&frames[f][k].r, &base);
                let perm: Option<Vec<usize>> = v.iter().map(|p| find_vertex(&v, mv(&r, *p))).collect();
                if let Some(p) = perm { syms.push(r); vperms.push(p); }
            }
        }
        let mut fperm = Vec::new();
        for p in &vperms {
            let mut fp = Vec::new();
            for f in &faces {
                let mut img: Vec<usize> = f.iter().map(|&i| p[i]).collect();
                img.sort();
                let j = faces.iter().position(|g| { let mut s = g.clone(); s.sort(); s == img }).unwrap();
                fp.push(j);
            }
            fperm.push(fp);
        }
        let axes = group_axes(&syms);
        let mut flags = Vec::new();
        for f in 0..faces.len() {
            for k in 0..faces[f].len() {
                let mut rep = true;
                for (s, p) in vperms.iter().enumerate() {
                    let f2 = fperm[s][f];
                    let k2 = faces[f2].iter().position(|&x| x == p[faces[f][k]]).unwrap();
                    if (f2, k2) < (f, k) { rep = false; break; }
                }
                if rep { flags.push((f, k)); }
            }
        }
        Solid { name, v, faces, fnorm, edges, edirs, radius, inr, syms, fperm, axes, frames, flags }
    }
}

pub fn load(path: &str) -> Vec<Solid> {
    let txt = std::fs::read_to_string(path).expect("solids file");
    let mut lines = txt.lines();
    let mut out = Vec::new();
    while let Some(h) = lines.next() {
        let p: Vec<&str> = h.split_whitespace().collect();
        if p.is_empty() { continue; }
        assert_eq!(p[0], "solid");
        let (nv, nf): (usize, usize) = (p[2].parse().unwrap(), p[3].parse().unwrap());
        let mut v = Vec::new();
        for _ in 0..nv {
            let c: Vec<f64> = lines.next().unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect();
            v.push([c[0], c[1], c[2]]);
        }
        let mut faces = Vec::new();
        for _ in 0..nf {
            let c: Vec<usize> = lines.next().unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect();
            faces.push(c[1..].to_vec());
        }
        out.push(Solid::new(p[1].to_string(), v, faces));
    }
    out
}
