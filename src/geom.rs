pub type V3 = [f64; 3];
pub type M3 = [[f64; 3]; 3];
pub const EPS: f64 = 1e-7;

pub fn add(a: V3, b: V3) -> V3 { [a[0] + b[0], a[1] + b[1], a[2] + b[2]] }
pub fn sub(a: V3, b: V3) -> V3 { [a[0] - b[0], a[1] - b[1], a[2] - b[2]] }
pub fn scale(a: V3, s: f64) -> V3 { [a[0] * s, a[1] * s, a[2] * s] }
pub fn dot(a: V3, b: V3) -> f64 { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }
pub fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
pub fn norm(a: V3) -> f64 { dot(a, a).sqrt() }
pub fn unit(a: V3) -> V3 { scale(a, 1.0 / norm(a)) }
pub fn dist(a: V3, b: V3) -> f64 { norm(sub(a, b)) }

pub const IDENT: M3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
pub fn mmul(a: &M3, b: &M3) -> M3 {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 { for j in 0..3 { r[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j]; } }
    r
}
pub fn mv(m: &M3, v: V3) -> V3 { [dot(m[0], v), dot(m[1], v), dot(m[2], v)] }
pub fn mt(m: &M3) -> M3 {
    [[m[0][0], m[1][0], m[2][0]], [m[0][1], m[1][1], m[2][1]], [m[0][2], m[1][2], m[2][2]]]
}
pub fn mclose(a: &M3, b: &M3, tol: f64) -> bool {
    for i in 0..3 { for j in 0..3 { if (a[i][j] - b[i][j]).abs() > tol { return false; } } }
    true
}
pub fn from_cols(a: V3, b: V3, c: V3) -> M3 {
    [[a[0], b[0], c[0]], [a[1], b[1], c[1]], [a[2], b[2], c[2]]]
}
/// Rotation by angle about unit axis u.
pub fn rot(u: V3, ang: f64) -> M3 { rot_sc(u, ang.sin_cos()) }
/// Rotation about unit axis u from the angle's (sin, cos).
pub fn rot_sc(u: V3, (s, c): (f64, f64)) -> M3 {
    let t = 1.0 - c;
    [
        [t * u[0] * u[0] + c, t * u[0] * u[1] - s * u[2], t * u[0] * u[2] + s * u[1]],
        [t * u[0] * u[1] + s * u[2], t * u[1] * u[1] + c, t * u[1] * u[2] - s * u[0]],
        [t * u[0] * u[2] - s * u[1], t * u[1] * u[2] + s * u[0], t * u[2] * u[2] + c],
    ]
}
/// Axis (unit, sign-normalised) of a non-identity rotation.
pub fn rot_axis(r: &M3) -> V3 {
    let a = [r[2][1] - r[1][2], r[0][2] - r[2][0], r[1][0] - r[0][1]];
    let mut u = if norm(a) > 1e-6 {
        unit(a)
    } else {
        // angle pi: R + I = 2 u u^T
        let mut best = [0.0; 3];
        let mut bn = 0.0;
        for j in 0..3 {
            let col = [r[0][j] + if j == 0 { 1.0 } else { 0.0 }, r[1][j] + if j == 1 { 1.0 } else { 0.0 }, r[2][j] + if j == 2 { 1.0 } else { 0.0 }];
            if norm(col) > bn { bn = norm(col); best = col; }
        }
        unit(best)
    };
    for i in 0..3 {
        if u[i].abs() > 1e-9 {
            if u[i] < 0.0 { u = scale(u, -1.0); }
            break;
        }
    }
    u
}

#[derive(Clone, Copy, Debug)]
pub struct Xf { pub r: M3, pub t: V3 }
impl Xf {
    pub fn ident() -> Xf { Xf { r: IDENT, t: [0.0; 3] } }
    pub fn ap(&self, p: V3) -> V3 { add(mv(&self.r, p), self.t) }
    pub fn mul(&self, o: &Xf) -> Xf { Xf { r: mmul(&self.r, &o.r), t: add(mv(&self.r, o.t), self.t) } }
    pub fn inv(&self) -> Xf { let rt = mt(&self.r); Xf { r: rt, t: scale(mv(&rt, self.t), -1.0) } }
    /// Rotation about the line through o with direction u.
    pub fn about(o: V3, u: V3, ang: f64) -> Xf { let r = rot(u, ang); Xf { r, t: sub(o, mv(&r, o)) } }
    pub fn about_sc(o: V3, u: V3, sc: (f64, f64)) -> Xf { let r = rot_sc(u, sc); Xf { r, t: sub(o, mv(&r, o)) } }
}

/// Close a set of generators into a finite group (elements as Xf); None if it exceeds cap.
pub fn closure(gens: &[Xf], cap: usize) -> Option<Vec<Xf>> {
    let mut el = vec![Xf::ident()];
    let mut i = 0;
    while i < el.len() {
        for g in gens {
            let n = g.mul(&el[i]);
            if !el.iter().any(|e| mclose(&e.r, &n.r, 1e-6)) {
                if el.len() >= cap { return None; }
                el.push(n);
            }
        }
        i += 1;
    }
    Some(el)
}

/// Distinct axes (dir, order) of a finite rotation group given by matrices.
pub fn group_axes(rs: &[M3]) -> Vec<(V3, usize)> {
    let mut axes: Vec<(V3, usize)> = Vec::new();
    for r in rs {
        if mclose(r, &IDENT, 1e-6) { continue; }
        let u = rot_axis(r);
        match axes.iter_mut().find(|(a, _)| dist(*a, u) < 1e-5) {
            Some(e) => e.1 += 1,
            None => axes.push((u, 2)),
        }
    }
    axes
}

/// Like closure, but returns the build recipe: element k = gens[a] * element[b].
pub fn closure_recipe(gens: &[Xf], cap: usize) -> Option<Vec<(usize, usize)>> {
    let mut el = vec![Xf::ident()];
    let mut rec = vec![(0, 0)];
    let mut i = 0;
    while i < el.len() {
        for (gi, g) in gens.iter().enumerate() {
            let n = g.mul(&el[i]);
            if !el.iter().any(|e| mclose(&e.r, &n.r, 1e-6)) {
                if el.len() >= cap { return None; }
                el.push(n);
                rec.push((gi, i));
            }
        }
        i += 1;
    }
    Some(rec)
}
/// Rebuild the group from a recipe into a reused buffer.
pub fn from_recipe(gens: &[Xf], rec: &[(usize, usize)], el: &mut Vec<Xf>) {
    el.clear();
    el.push(Xf::ident());
    for &(g, b) in &rec[1..] { let n = gens[g].mul(&el[b]); el.push(n); }
}
