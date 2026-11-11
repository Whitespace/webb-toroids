//! Data-parallel kernels over a solid's vertices in structure-of-arrays form.
//!
//! No intrinsics: every loop works on fixed blocks of `W` lanes with independent per-lane
//! state, which LLVM turns into NEON on aarch64 and AVX2/AVX-512 on x86-64 (build with
//! `-C target-cpu=native`, see .cargo/config.toml). The arrays are padded to a multiple of `W`
//! by repeating vertex 0, so there are no remainder loops and padding never changes a
//! max/min/any result.
use crate::geom::V3;

pub const W: usize = 8;

pub struct Soa {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    pub n: usize,
}

impl Soa {
    pub fn new(v: &[V3]) -> Soa {
        let n = v.len();
        let m = (n + W - 1) / W * W;
        let col = |i: usize| -> Vec<f64> { (0..m).map(|k| if k < n { v[k][i] } else { v[0][i] }).collect() };
        Soa { x: col(0), y: col(1), z: col(2), n }
    }

    /// max over vertices of a . v
    #[inline]
    pub fn proj_max(&self, a: V3) -> f64 {
        let mut m = [f64::MIN; W];
        for ((x, y), z) in self.x.chunks_exact(W).zip(self.y.chunks_exact(W)).zip(self.z.chunks_exact(W)) {
            for l in 0..W {
                let d = a[0] * x[l] + a[1] * y[l] + a[2] * z[l];
                m[l] = if d > m[l] { d } else { m[l] };
            }
        }
        let mut r = m[0];
        for l in 1..W { if m[l] > r { r = m[l]; } }
        r
    }

    /// max over vertices of |v - o|^2
    #[inline]
    pub fn far2(&self, o: V3) -> f64 {
        let mut m = [f64::MIN; W];
        for ((x, y), z) in self.x.chunks_exact(W).zip(self.y.chunks_exact(W)).zip(self.z.chunks_exact(W)) {
            for l in 0..W {
                let d = d2(o, x[l], y[l], z[l]);
                m[l] = if d > m[l] { d } else { m[l] };
            }
        }
        let mut r = m[0];
        for l in 1..W { if m[l] > r { r = m[l]; } }
        r
    }

    /// true if some vertex is at distance 1 (to tol on the square) from o
    #[inline]
    pub fn any_unit(&self, o: V3, tol: f64) -> bool {
        let mut any = false;
        for ((x, y), z) in self.x.chunks_exact(W).zip(self.y.chunks_exact(W)).zip(self.z.chunks_exact(W)) {
            let mut hit = [false; W];
            for l in 0..W { hit[l] = (d2(o, x[l], y[l], z[l]) - 1.0).abs() < tol; }
            for l in 0..W { any |= hit[l]; }
        }
        any
    }

    /// Wedge test for the supporting-plane scan. Each vertex projects to q = (a.v + a0, b.v + b0).
    /// Returns false only if no vertex could change the scan's state: every q lies inside the
    /// wedge from l round to h (cross(h,q) <= 0, cross(l,q) >= 0) and none is antiparallel to r.
    /// A side that has no bound yet (`has_h` / `has_l` false) must have no vertex at all.
    #[inline]
    pub fn wedge_out(&self, a: V3, a0: f64, b: V3, b0: f64, h: (f64, f64), has_h: bool, lo: (f64, f64), has_l: bool, r: (f64, f64)) -> bool {
        let k = 1e-13 * (r.0 * r.0 + r.1 * r.1);
        for ((x, y), z) in self.x.chunks_exact(W).zip(self.y.chunks_exact(W)).zip(self.z.chunks_exact(W)) {
            let mut bad = [false; W];
            for l in 0..W {
                let q0 = proj(a, a0, x[l], y[l], z[l]);
                let q1 = proj(b, b0, x[l], y[l], z[l]);
                let ch = h.0 * q1 - h.1 * q0;
                let cl = lo.0 * q1 - lo.1 * q0;
                let s = r.0 * q1 - r.1 * q0;
                let col = s * s <= k * (q0 * q0 + q1 * q1);
                let opp = col & (r.0 * q0 + r.1 * q1 < 0.0);
                let bh = if has_h { ch > 0.0 } else { !col & (s > 0.0) };
                let bl = if has_l { cl < 0.0 } else { !col & (s < 0.0) };
                bad[l] = bh | bl | opp;
            }
            let mut any = false;
            for l in 0..W { any |= bad[l]; }
            if any { return true; }
        }
        false
    }
}

#[inline(always)]
pub fn proj(a: V3, a0: f64, x: f64, y: f64, z: f64) -> f64 { a[0] * x + a[1] * y + a[2] * z + a0 }

#[inline(always)]
pub fn d2(o: V3, x: f64, y: f64, z: f64) -> f64 {
    let (dx, dy, dz) = (x - o[0], y - o[1], z - o[2]);
    dx * dx + dy * dy + dz * dz
}

/// out[i] = u . w[i] over SoA direction tables (lengths padded to a multiple of W).
#[inline]
pub fn dots(u: V3, wx: &[f64], wy: &[f64], wz: &[f64], out: &mut [f64]) {
    for (((x, y), z), o) in wx.chunks_exact(W).zip(wy.chunks_exact(W)).zip(wz.chunks_exact(W)).zip(out.chunks_exact_mut(W)) {
        for l in 0..W { o[l] = u[0] * x[l] + u[1] * y[l] + u[2] * z[l]; }
    }
}
