mod geom;
mod simd;
mod solids;
mod validate;
use geom::*;
use solids::Solid;
use std::collections::HashSet;
use std::io::Write;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;

static G_TESTS: AtomicU64 = AtomicU64::new(0);
static G_HITS: AtomicU64 = AtomicU64::new(0);
static G_QUICK: AtomicU64 = AtomicU64::new(0);
static G_VAL_NS: AtomicU64 = AtomicU64::new(0);
const RNAMES: [&str; 12] = ["overlap-centre", "overlap-chain", "overlap-sphere", "overlap", "contact", "q-prescreen", "q-proxy", "nonmanifold", "pinch", "disconnected", "genus0", "other"];
static G_R: [AtomicU64; 12] = [const { AtomicU64::new(0) }; 12];
static G_RNS: [AtomicU64; 12] = [const { AtomicU64::new(0) }; 12];
static G_DONE: AtomicU64 = AtomicU64::new(0);

static LOG: Mutex<Option<std::fs::File>> = Mutex::new(None);
/// Progress line to stderr and, once the output directory is open, to OUT/run.log.
fn say(msg: String) {
    eprintln!("{}", msg);
    if let Some(f) = LOG.lock().unwrap().as_mut() { let _ = writeln!(f, "{}", msg); }
}

const FLIP: M3 = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]];

struct AxEntry { solid: usize, face: usize, k: usize, order: usize, xp: Xf, p: V3, w: V3, bq: [usize; 3], nbq: usize }
#[derive(Clone)]
struct Link { solid: usize, g: usize, k: usize, h: usize, hn: usize, xp: Xf, exit: Xf, through: bool }

/// `axw[n]` holds the axis directions of `ax[n]` as padded structure-of-arrays for the vector dot product.
/// Distinct axis directions of `ax[n]` (padded structure-of-arrays) and, per direction, the entries that have it.
#[derive(Default)]
struct UDirs { w: [Vec<f64>; 3], n: usize, start: Vec<u32>, ent: Vec<u32>, p: Vec<V3> }
/// Links of one list that share an exit rotation (and exit polygon): the seed axis direction in
/// their exit frames is the same, so one scan of the direction table serves all of them.
struct RotClass { r: M3, hn: usize, links: Vec<u32>, t: Vec<V3> }

/// `pqk[i]` goes with `pq[i]`: a top-order axis of the generated group nearest the first generator's
/// axis, in the frame (e2, e3, u) of the two generator axes, and how many such axes ring it.
struct Tables { pqk: Vec<(V3, usize)>, ud: Vec<UDirs>, rc: [Vec<Vec<RotClass>>; 2], axw: Vec<[Vec<f64>; 3]>, sc: [(f64, f64); 6], ax: Vec<Vec<AxEntry>>, links: Vec<Vec<Link>>, thru: Vec<Vec<Link>>, cosset: Vec<f64>, pq: Vec<(usize, usize, f64, Vec<(usize, usize)>)>, cosbin: Vec<bool> }

fn build_tables(lib: &[Solid], active: &[bool], linkok: &[bool], struts: bool) -> Tables {
    let mut ax: Vec<Vec<AxEntry>> = (0..11).map(|_| Vec::new()).collect();
    let mut links: Vec<Vec<Link>> = (0..11).map(|_| Vec::new()).collect();
    for (si, s) in lib.iter().enumerate() {
        if !active[si] { continue; }
        for &(g, k) in &s.flags {
            let n = s.faces[g].len();
            let fi = s.frames[g][k].inv();
            let xp = Xf { r: mmul(&FLIP, &fi.r), t: mv(&FLIP, fi.t) };
            for &(u, order) in &s.axes {
                // orders usable for B's rotation: the divisors >= 3, else the 2-fold
                let mut bq = [2usize; 3];
                let mut nbq = 0;
                for d in divisors(order) { if d >= 3 { bq[nbq] = d; nbq += 1; } }
                if nbq == 0 { nbq = 1; }
                ax[n].push(AxEntry { solid: si, face: g, k, order, xp, p: xp.t, w: mv(&xp.r, u), bq, nbq });
            }
            for h in 0..s.faces.len() {
                if h == g || !linkok[si] { continue; }
                let through = dot(s.fnorm[h], s.fnorm[g]) < -1.0 + 1e-9;
                // strut mode: the free ("through") links are only base to base through a prism or
                // antiprism over this polygon; every other link counts against --any
                let through = if struts { through && s.faces[h].len() == n && s.v.len() == 2 * n } else { through };
                links[n].push(Link { solid: si, g, k, h, hn: s.faces[h].len(), xp, exit: xp.mul(&s.frames[h][0]), through });
            }
        }
    }
    // allowed |cos| between distinct axes of the polyhedral groups
    let phi = (1.0 + 5f64.sqrt()) / 2.0;
    let o0 = [0.0; 3];
    let gi = closure(&[Xf::about(o0, unit([0.0, 1.0, phi]), 0.4 * std::f64::consts::PI), Xf::about(o0, unit([1.0, 1.0, 1.0]), 2.0 * std::f64::consts::PI / 3.0)], 60).unwrap();
    let go = closure(&[Xf::about(o0, [0.0, 0.0, 1.0], 0.5 * std::f64::consts::PI), Xf::about(o0, [1.0, 0.0, 0.0], 0.5 * std::f64::consts::PI)], 24).unwrap();
    assert_eq!(gi.len(), 60);
    assert_eq!(go.len(), 24);
    let mut cosset: Vec<f64> = Vec::new();
    let mut pq: Vec<(usize, usize, f64, Vec<(usize, usize)>)> = Vec::new();
    for g in [&gi, &go] {
        let rs: Vec<M3> = g.iter().map(|x| x.r).collect();
        let axs = group_axes(&rs);
        for i in 0..axs.len() { for j in 0..axs.len() { for sg in [1.0, -1.0] {
            if i == j { continue; }
            let aj = scale(axs[j].0, sg);
            let c = dot(axs[i].0, aj);
            if !cosset.iter().any(|x| (x - c.abs()).abs() < 1e-9) { cosset.push(c.abs()); }
            if pq.iter().any(|x| x.0 == axs[i].1 && x.1 == axs[j].1 && (x.2 - c).abs() < 1e-9) { continue; }
            // signed cos fixes the pair of generators up to a rotation, so the recipe transfers
            let o0 = [0.0; 3];
            let tau = 2.0 * std::f64::consts::PI;
            if let Some(rec) = closure_recipe(&[Xf::about(o0, axs[i].0, tau / axs[i].1 as f64), Xf::about(o0, aj, tau / axs[j].1 as f64)], 60) {
                pq.push((axs[i].1, axs[j].1, c, rec));
            }
        } } }
    }
    let mut cosbin = vec![false; 4098];
    for c in &cosset { let b = (c * 4096.0) as usize; for k in b.saturating_sub(1)..=(b + 1).min(4097) { cosbin[k] = true; } }
    let thru: Vec<Vec<Link>> = links.iter().map(|v| v.iter().filter(|l| l.through).cloned().collect()).collect();
    let axw: Vec<[Vec<f64>; 3]> = ax.iter().map(|v| {
        let m = (v.len() + simd::W - 1) / simd::W * simd::W;
        let col = |i: usize| -> Vec<f64> { (0..m).map(|k| if k < v.len() { v[k].w[i] } else { 0.0 }).collect() };
        [col(0), col(1), col(2)]
    }).collect();
    let mut sc = [(0.0, 1.0); 6];
    for k in 2..6 { sc[k] = (2.0 * std::f64::consts::PI / k as f64).sin_cos(); }
    let ud: Vec<UDirs> = ax.iter().map(|v| {
        let mut keys: std::collections::HashMap<[i64; 3], usize> = Default::default();
        let mut dirs: Vec<V3> = Vec::new();
        let mut mem: Vec<Vec<u32>> = Vec::new();
        for (i, e) in v.iter().enumerate() {
            let k = [(e.w[0] * 1e8).round() as i64, (e.w[1] * 1e8).round() as i64, (e.w[2] * 1e8).round() as i64];
            let d = *keys.entry(k).or_insert_with(|| { dirs.push(e.w); mem.push(Vec::new()); dirs.len() - 1 });
            mem[d].push(i as u32);
        }
        let m = (dirs.len() + simd::W - 1) / simd::W * simd::W;
        let col = |i: usize| -> Vec<f64> { (0..m).map(|k| if k < dirs.len() { dirs[k][i] } else { 0.0 }).collect() };
        let mut start = vec![0u32];
        let mut ent = Vec::new();
        for g in &mem { ent.extend_from_slice(g); start.push(ent.len() as u32); }
        let p = ent.iter().map(|&i| v[i as usize].p).collect();
        UDirs { w: [col(0), col(1), col(2)], n: dirs.len(), start, ent, p }
    }).collect();
    let classes = |lists: &Vec<Vec<Link>>| -> Vec<Vec<RotClass>> { lists.iter().map(|v| {
        let mut keys: std::collections::HashMap<[i64; 10], usize> = Default::default();
        let mut out: Vec<RotClass> = Vec::new();
        for (i, l) in v.iter().enumerate() {
            let mut k = [l.hn as i64; 10];
            for a in 0..3 { for b in 0..3 { k[a * 3 + b] = (l.exit.r[a][b] * 1e8).round() as i64; } }
            let c = *keys.entry(k).or_insert_with(|| { out.push(RotClass { r: l.exit.r, hn: l.hn, links: Vec::new(), t: Vec::new() }); out.len() - 1 });
            out[c].links.push(i as u32);
            out[c].t.push(l.exit.t);
        }
        out
    }).collect() };
    let rc = [classes(&thru), classes(&links)];
    let tau = 2.0 * std::f64::consts::PI;
    let pqk: Vec<(V3, usize)> = pq.iter().map(|(p, q, c, rec)| {
        let gens = [Xf::about(o0, [0.0, 0.0, 1.0], tau / *p as f64), Xf::about(o0, [(1.0 - c * c).sqrt(), 0.0, *c], tau / *q as f64)];
        let mut el = Vec::new();
        from_recipe(&gens, rec, &mut el);
        let rs: Vec<M3> = el.iter().map(|x| x.r).collect();
        let mut dirs: Vec<(V3, usize)> = Vec::new();
        for (d, o) in group_axes(&rs) { if d[2].abs() < 1.0 - 1e-6 { dirs.push((d, o)); dirs.push((scale(d, -1.0), o)); } }
        let top = dirs.iter().map(|x| x.1).max().unwrap_or(0);
        let cmax = dirs.iter().filter(|x| x.1 == top).map(|x| x.0[2]).fold(-2.0, f64::max);
        let sel: Vec<V3> = dirs.iter().filter(|x| x.1 == top && x.0[2] > cmax - 1e-6).map(|x| x.0).collect();
        let m = sel.len();
        if m == 0 { return ([0.0; 3], 0); }
        // they must be evenly spaced about the first axis, or the key below is not canonical
        let mut az: Vec<f64> = sel.iter().map(|d| d[1].atan2(d[0]).rem_euclid(tau / m as f64)).collect();
        az.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let spread = az[m - 1] - az[0];
        if spread > 1e-6 && (tau / m as f64 - spread) > 1e-6 { return ([0.0; 3], 0); }
        (sel[0], m)
    }).collect();
    Tables { pqk, ud, rc, axw, sc, ax, links, thru, cosset, pq, cosbin }
}

struct Shared {
    seen: HashSet<String>,
    tsv: std::fs::File,
    jsonl: std::fs::File,
    ckpt: std::fs::File,
    n: usize,
    reasons: std::collections::HashMap<&'static str, u64>,
}

struct Ctx<'a> {
    lib: &'a [Solid],
    tab: &'a Tables,
    a: usize,
    u: V3,
    maxdepth: usize,
    ordered: bool,
    minorder: usize,
    near: usize,
    shared: &'a Mutex<Shared>,
    outdir: &'a str,
    tests: u64,
    hits: u64,
    desc: Vec<String>,
    reasons: std::collections::HashMap<&'static str, u64>,
    scratch: validate::Scratch,
    /// divisors >= 3 of the seed axis order
    pdiv: Vec<usize>,
    /// a centre s*u with s strictly inside this interval is strictly inside the seed piece
    inside: (f64, f64),
    csbuf: Vec<f64>,
    st: Stats,
    /// rotations about the seed axis, one per usable order; known before the chain is built
    axrot: Vec<M3>,
    /// --max-radius: no piece centre may be further than this from the seed axis or the centre
    maxrad: f64,
    pruned: u64,
    /// join the last chain level through the direction tables instead of scanning link by link
    join: bool,
    /// bitmap over the signed cos of the angle between the seed axis and an end axis: the values
    /// some group recipe accepts for this seed's axis orders, one bin of slack either side
    sbin: Box<[bool; 65536]>,
    bins: Vec<i32>,
    /// the signed cos values behind `sbin`
    scos: Vec<f64>,
    /// Per-node cache: group (as a set) -> what the chain so far, without link or end piece
    /// still to come, does under it. Cleared in `test`, which every node runs first.
    share: bool,
    check: bool,
    gc: std::collections::HashMap<[i64; 3], u32, std::hash::BuildHasherDefault<validate::Fx>>,
    gce: Vec<GEnt>,
    gcg: Vec<Xf>,
    gc_np: usize,
    ex: V3,
    ey: V3,
    pre_ov: u64,
    pre_q: u64,
    cand: Vec<(u32, u32)>,
}

/// Per-thread counters, pushed to the shared atomics now and then so the hot path never
/// touches a cache line another thread writes.
#[derive(Default)]
struct Stats { calls: u64, tests: u64, hits: u64, quick: u64, quick_total: u64, val_ns: u64, r: [u64; 12], rns: [u64; 12] }
impl Stats {
    fn flush(&mut self) {
        G_TESTS.fetch_add(self.tests, Ordering::Relaxed);
        G_HITS.fetch_add(self.hits, Ordering::Relaxed);
        G_QUICK.fetch_add(self.quick, Ordering::Relaxed);
        G_VAL_NS.fetch_add(self.val_ns, Ordering::Relaxed);
        for k in 0..12 { if self.r[k] > 0 { G_R[k].fetch_add(self.r[k], Ordering::Relaxed); G_RNS[k].fetch_add(self.rns[k], Ordering::Relaxed); } }
        self.quick_total += self.quick;
        *self = Stats { calls: self.calls, quick_total: self.quick_total, ..Default::default() };
    }
}

/// verdict 1: the prefix orbit overlaps itself or is too big; 2: it fails the outermost-vertex
/// screen at `vs`; 0: neither.
struct GEnt { verdict: u8, vs: V3, far: f64, g0: usize, gl: usize }

fn divisors(n: usize) -> Vec<usize> { [2usize, 3, 4, 5].iter().cloned().filter(|d| n % d == 0).collect() }

fn xf_json(x: &Xf) -> String {
    format!("[{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}]",
        x.r[0][0], x.r[0][1], x.r[0][2], x.r[1][0], x.r[1][1], x.r[1][2], x.r[2][0], x.r[2][1], x.r[2][2], x.t[0], x.t[1], x.t[2])
}

pub fn write_off(path: &str, r: &validate::Res, lib: &[Solid]) {
    // companion file: the solid each boundary face belongs to
    let pcs: Vec<&str> = r.fsolid.iter().map(|s| lib[*s].name.as_str()).collect();
    std::fs::write(format!("{}.pieces", path), pcs.join("\n")).unwrap();
    let mut s = format!("OFF\n{} {} 0\n", r.pts.len(), r.faces.len());
    for p in &r.pts { s += &format!("{:.9} {:.9} {:.9}\n", p[0], p[1], p[2]); }
    for f in &r.faces {
        s += &format!("{}", f.len());
        for i in f { s += &format!(" {}", i); }
        s += "\n";
    }
    std::fs::write(path, s).unwrap();
}

fn record(ctx: &mut Ctx, reps: &[(usize, Xf)], gens: &[Xf], centre: V3, res: Result<validate::Res, &'static str>) {
    let r = match res {
        Err(e) => { *ctx.reasons.entry(e).or_insert(0) += 1; return; }
        Ok(r) => r,
    };
    *ctx.reasons.entry("valid").or_insert(0) += 1;
    let good = r.a_viol == 0 && r.q_ok;
    if good { *ctx.reasons.entry("RAQT").or_insert(0) += 1; }
    if std::env::var("COUNT_UNIQUE").is_ok() {
        let (mut m1, mut m3) = (0.0, 0.0);
        for p in &r.pts { let d = dist(*p, centre); m1 += d; m3 += d * d * d; }
        let sig = format!("{} {} {} {} {} {} {} {:.4} {:.3} {:.2}", r.gorder, r.genus, r.nv, r.ne, r.nf, r.a_viol, r.q_missing, r.rcirc, m1, m3);
        if ctx.shared.lock().unwrap().seen.insert(sig) { *ctx.reasons.entry("valid-unique").or_insert(0) += 1; }
    }
    if !(good || (r.a_viol == 0 && !r.hull_faces.is_empty() && r.q_missing <= ctx.near)) { return; }
    let mut sh = ctx.shared.lock().unwrap();
    let names: Vec<String> = r.counts.iter().map(|(s, c)| format!("{}x{}", c, ctx.lib[*s].name)).collect();
    let (mut m1, mut m3) = (0.0, 0.0);
    for p in &r.pts { let d = dist(*p, centre); m1 += d; m3 += d * d * d; }
    let sig = format!("{} {} {} {} {} {} {} {:.4} {:.3} {:.2} {}", r.gorder, r.genus, r.nv, r.ne, r.nf, r.a_viol, r.q_missing, r.rcirc, m1, m3, names.join("+"));
    if !sh.seen.insert(sig.clone()) { return; }
    let id = sh.n;
    sh.n += 1;

    let hf: Vec<String> = r.hull_faces.iter().map(|n| n.to_string()).collect();
    writeln!(sh.tsv, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}",
        id, if good { "RAQT" } else { "-" }, r.gorder, r.genus, r.npieces, r.a_viol, r.q_missing, r.q_edges,
        if r.hull_regular { "regular" } else { "irregular" }, r.rin, r.rcirc, if r.rin < 1e9 { r.rin / r.rcirc } else { 0.0 },
        hf.join(","), r.nv, r.ne, r.nf, names.join("+"), ctx.desc.join(" > ")).unwrap();
    let rj: Vec<String> = reps.iter().map(|(s, x)| format!("{{\"solid\":\"{}\",\"xf\":{}}}", ctx.lib[*s].name, xf_json(x))).collect();
    let gj: Vec<String> = gens.iter().map(xf_json).collect();
    writeln!(sh.jsonl, "{{\"id\":{},\"good\":{},\"group\":{},\"genus\":{},\"centre\":[{:.15e},{:.15e},{:.15e}],\"reps\":[{}],\"gens\":[{}]}}",
        id, good, r.gorder, r.genus, centre[0], centre[1], centre[2], rj.join(","), gj.join(",")).unwrap();
    write_off(&format!("{}/off/{:06}.off", ctx.outdir, id), &r, ctx.lib);
    // last, so a result is only treated as known on resume once all its files exist
    writeln!(sh.ckpt, "sig {}\t{}", id, sig).unwrap();
}

fn test(ctx: &mut Ctx, w: &Xf, n: usize, pieces: &mut Vec<(usize, Xf)>) {
    if !ctx.gce.is_empty() { ctx.gc.clear(); ctx.gce.clear(); ctx.gcg.clear(); }
    ctx.gc_np = pieces.len();
    let wi = w.inv();
    let o = wi.t;
    let u = mv(&wi.r, ctx.u);
    let tab = ctx.tab;
    let ents = &tab.ax[n];
    ctx.tests += ents.len() as u64;
    ctx.st.tests += ents.len() as u64;
    ctx.st.calls += 1;
    if ctx.st.calls & 0xfff == 0 { ctx.st.flush(); }
    // cos of the angle between the seed axis and every candidate axis, in one vector pass
    let mut buf = std::mem::take(&mut ctx.csbuf);
    let [wx, wy, wz] = &tab.axw[n];
    buf.resize(wx.len(), 0.0);
    simd::dots(u, wx, wy, wz, &mut buf);
    for i in 0..ents.len() {
        let cs = buf[i];
        let c = cs.abs();
        if !tab.cosbin[(c * 4096.0) as usize] { continue; }
        entry(ctx, w, o, u, n, i, cs, pieces);
    }
    ctx.csbuf = buf;
}

/// One end piece: entry `i` of `ax[n]` on the open face with frame `w`. `o`, `u` are the seed axis
/// (a point and the direction) in that frame, `cs` the cos of the angle between the two axes.
fn entry(ctx: &mut Ctx, w: &Xf, o: V3, u: V3, n: usize, i: usize, cs: f64, pieces: &mut Vec<(usize, Xf)>) {
    let tab = ctx.tab;
    let c = cs.abs();
    {
        if !tab.cosset.iter().any(|x| (x - c).abs() < 1e-7) { return; }
        let e = &tab.ax[n][i];
        let nn = cross(u, e.w);
        let d = sub(e.p, o);
        if dot(d, nn).abs() > 1e-7 { return; }
        let bq = &e.bq[..e.nbq];
        if ctx.ordered && bq[0] >= 3 && e.solid < ctx.a { return; }
        let s = dot(cross(d, e.w), nn) / dot(nn, nn);
        let centre = scale(ctx.u, s);
        // a piece sitting on the centre gives a star of arms, not a shell; skip those
        let db = dist(centre, w.ap(e.p));
        if s.abs() < 1e-6 || db < 1e-6 { return; }
        if s.abs() > ctx.maxrad || db > ctx.maxrad { return; }
        // centre strictly inside the seed piece: it would meet its own images
        if s > ctx.inside.0 && s < ctx.inside.1 { return; }
        let xb = w.mul(&e.xp);
        let wb = mv(&w.r, e.w);
        pieces.push((e.solid, xb));
        for pi in 0..ctx.pdiv.len() {
            let p = ctx.pdiv[pi];
            for &q in bq {
                let ri = match tab.pq.iter().position(|x| x.0 == p && x.1 == q && (x.2 - cs).abs() < 1e-7 && x.3.len() >= ctx.minorder) { Some(x) => x, None => continue };
                let rec = &tab.pq[ri].3;
                ctx.hits += 1;
                ctx.st.hits += 1;
                let gens = [Xf::about_sc(centre, ctx.u, tab.sc[p]), Xf::about_sc(centre, wb, tab.sc[q])];
                // Many hits at one node generate the same group, and most die on what the chain
                // prefix does under it, whatever comes after. So that is worked out once per group.
                let (rv, m) = tab.pqk[ri];
                if ctx.share && m > 0 {
                    let e2 = unit(sub(wb, scale(ctx.u, cs)));
                    let e3 = cross(ctx.u, e2);
                    let r = add(add(scale(e2, rv[0]), scale(e3, rv[1])), scale(ctx.u, rv[2]));
                    let per = 2.0 * std::f64::consts::PI / m as f64;
                    let az = dot(r, ctx.ey).atan2(dot(r, ctx.ex)).rem_euclid(per);
                    let key = [(rec.len() * 8 + p) as i64, (s * 1e6).round() as i64, (az * 1e6).round() as i64];
                    let gi = match ctx.gc.get(&key) {
                        Some(&i) => i as usize,
                        None => {
                            from_recipe(&gens, rec, &mut ctx.scratch.grp);
                            ctx.scratch.probe = true;
                            let res = validate::validate(ctx.lib, &pieces[..ctx.gc_np], centre, true, &mut ctx.scratch);
                            ctx.scratch.probe = false;
                            let verdict = match res { Err("overlap-centre") | Err("overlap-sphere") | Err("too-big") => 1, Err("q-prescreen") => 2, _ => 0 };
                            let g0 = ctx.gcg.len();
                            if verdict == 2 || ctx.check { ctx.gcg.extend_from_slice(&ctx.scratch.grp); }
                            ctx.gce.push(GEnt { verdict, vs: ctx.scratch.vs, far: ctx.scratch.far, g0, gl: ctx.gcg.len() - g0 });
                            ctx.gc.insert(key, ctx.gce.len() as u32 - 1);
                            ctx.gce.len() - 1
                        }
                    };
                    let ge = &ctx.gce[gi];
                    let grp = &ctx.gcg[ge.g0..ge.g0 + ge.gl];
                    if ctx.check {
                        let mut el = Vec::new();
                        from_recipe(&gens, rec, &mut el);
                        assert!(el.len() == grp.len() && el.iter().all(|a| grp.iter().any(|b| mclose(&a.r, &b.r, 1e-6) && dist(a.t, b.t) < 1e-6)), "group cache key is not canonical");
                    }
                    if ge.verdict == 1 { ctx.pre_ov += 1; continue; }
                    if ge.verdict == 2 && ctx.near == 0 && pieces[ctx.gc_np..].iter().all(|(s, x)| validate::extra_clear(ctx.lib, *s, x, centre, grp, ge.vs, ge.far)) { ctx.pre_q += 1; continue; }
                }
                // quick reject: inscribed spheres of chain pieces against their images under the generators
                let mut quick = false;
                'q: for g in &gens {
                    for &(s1, x1) in pieces.iter() {
                        let c = g.ap(x1.t);
                        for &(s2, x2) in pieces.iter() {
                            let r = ctx.lib[s1].inr + ctx.lib[s2].inr - 1e-4;
                            let d = sub(c, x2.t);
                            let d2 = dot(d, d);
                            if d2 > 1e-10 && d2 < r * r { quick = true; break 'q; }
                        }
                    }
                }
                if quick { ctx.st.quick += 1; continue; }
                let t0 = std::time::Instant::now();
                from_recipe(&gens, rec, &mut ctx.scratch.grp);
                let res = validate::validate(ctx.lib, pieces, centre, ctx.near == 0, &mut ctx.scratch);
                let ns = t0.elapsed().as_nanos() as u64;
                ctx.st.val_ns += ns;
                if let Err(r) = &res { let k = RNAMES.iter().position(|x| x == r).unwrap_or(11); ctx.st.r[k] += 1; ctx.st.rns[k] += ns; }
                if res.is_ok() {
                    ctx.desc.push(format!("{}[f{} k{} ax{}]", ctx.lib[e.solid].name, e.face, e.k, e.order));
                    let reps: Vec<(usize, Xf)> = pieces.clone();
                    record(ctx, &reps, &gens, centre, res);
                    ctx.desc.pop();
                } else if let Err(r) = res { *ctx.reasons.entry(r).or_insert(0) += 1; }
            }
        }
        pieces.pop();
    }
}

/// True if no chain through this new link can survive, so the whole subtree is skipped. All
/// three tests are ones every candidate below would fail later anyway:
/// the link is too far from the seed axis (--max-radius); its inscribed sphere overlaps that of
/// an earlier chain piece; or, for every usable rotation about the seed axis, it overlaps an
/// image of a chain piece (the same test as the quick reject, done before the subtree is built).
fn prune(ctx: &Ctx, sl: usize, cl: V3, pieces: &[(usize, Xf)]) -> bool {
    let u = ctx.u;
    let off = sub(cl, scale(u, dot(cl, u)));
    if dot(off, off) > ctx.maxrad * ctx.maxrad { return true; }
    let il = ctx.lib[sl].inr;
    for &(s, x) in pieces {
        let r = il + ctx.lib[s].inr - 1e-4;
        let d = sub(cl, x.t);
        let d2 = dot(d, d);
        if d2 > 1e-10 && d2 < r * r { return true; }
    }
    if ctx.axrot.is_empty() { return false; }
    'rot: for g in &ctx.axrot {
        // the link against its own image
        let c = mv(g, cl);
        let d = sub(c, cl);
        let d2 = dot(d, d);
        let r = 2.0 * il - 1e-4;
        if d2 > 1e-10 && d2 < r * r { continue 'rot; }
        for &(s, x) in pieces {
            let r = il + ctx.lib[s].inr - 1e-4;
            // image of the link against the piece, and image of the piece against the link
            let d = sub(c, x.t);
            let d2 = dot(d, d);
            if d2 > 1e-10 && d2 < r * r { continue 'rot; }
            let d = sub(mv(g, x.t), cl);
            let d2 = dot(d, d);
            if d2 > 1e-10 && d2 < r * r { continue 'rot; }
        }
        return false;
    }
    true
}

/// The last chain level, without visiting every link. Gives exactly what one more level of
/// `dfs` would: the same (link, end piece) pairs reach `entry` in the same order. Links that
/// share an exit rotation see the seed axis in the same direction, so the direction table of the
/// exit polygon is scanned once per rotation class; only the few pairs whose axes make an
/// allowed angle and are coplanar (loose tolerances here, the exact ones in `entry`) go further.
fn join(ctx: &mut Ctx, w: &Xf, n: usize, any_left: usize, pieces: &mut Vec<(usize, Xf)>) {
    let tab = ctx.tab;
    let (list, classes) = if any_left == 0 { (&tab.thru[n], &tab.rc[0][n]) } else { (&tab.links[n], &tab.rc[1][n]) };
    let wi = w.inv();
    let o = wi.t;
    let u = mv(&wi.r, ctx.u);
    let mut bins = std::mem::take(&mut ctx.bins);
    let sbin: &[bool; 65536] = &ctx.sbin;
    let scos = std::mem::take(&mut ctx.scos);
    let mut cand = std::mem::take(&mut ctx.cand);
    cand.clear();
    let mut work = 0u64;
    for rc in classes {
        // the seed axis direction in the exit frame: transpose(r) * u
        let ul = [u[0] * rc.r[0][0] + u[1] * rc.r[1][0] + u[2] * rc.r[2][0], u[0] * rc.r[0][1] + u[1] * rc.r[1][1] + u[2] * rc.r[2][1], u[0] * rc.r[0][2] + u[1] * rc.r[1][2] + u[2] * rc.r[2][2]];
        let ud = &tab.ud[rc.hn];
        let [wx, wy, wz] = &ud.w;
        bins.resize(wx.len(), 0);
        simd::bins(ul, wx, wy, wz, &mut bins);
        work += ud.n as u64;
        for (ci, ch) in bins.chunks_exact(simd::W).enumerate() {
          // one branch per chunk: almost every chunk has no direction at an allowed angle
          if !ch.iter().fold(false, |a, &b| a | sbin[b as usize & 0xffff]) { continue; }
          for d in ci * simd::W..((ci + 1) * simd::W).min(ud.n) {
            if !sbin[bins[d] as usize & 0xffff] { continue; }
            // the exact angle test (the bins are coarse), then coplanarity on compact copies
            let wd = [wx[d], wy[d], wz[d]];
            let cs = dot(ul, wd);
            if !scos.iter().any(|c| (c - cs).abs() < 1e-6) { continue; }
            let nn = cross(ul, wd);
            let (m0, m1) = (ud.start[d] as usize, ud.start[d + 1] as usize);
            for (k, &li) in rc.links.iter().enumerate() {
                // the seed axis point in this link's exit frame, along nn
                let t = sub(o, rc.t[k]);
                let ol = [t[0] * rc.r[0][0] + t[1] * rc.r[1][0] + t[2] * rc.r[2][0], t[0] * rc.r[0][1] + t[1] * rc.r[1][1] + t[2] * rc.r[2][1], t[0] * rc.r[0][2] + t[1] * rc.r[1][2] + t[2] * rc.r[2][2]];
                let h = dot(ol, nn);
                for m in m0..m1 {
                    if (dot(ud.p[m], nn) - h).abs() < 1e-5 { cand.push((li, ud.ent[m])); }
                }
            }
          }
        }
    }
    ctx.tests += work;
    ctx.st.tests += work;
    ctx.st.calls += 1;
    if ctx.st.calls & 0xff == 0 { ctx.st.flush(); }
    cand.sort_unstable();
    let mut i = 0;
    while i < cand.len() {
        let li = cand[i].0;
        let mut j = i;
        while j < cand.len() && cand[j].0 == li { j += 1; }
        let l = &list[li as usize];
        let xl = w.mul(&l.xp);
        if prune(ctx, l.solid, xl.t, pieces) { ctx.pruned += 1; i = j; continue; }
        pieces.push((l.solid, xl));
        ctx.desc.push(format!("{}[f{} k{} > f{}]", ctx.lib[l.solid].name, l.g, l.k, l.h));
        let w2 = w.mul(&l.exit);
        let wi2 = w2.inv();
        let (o2, u2) = (wi2.t, mv(&wi2.r, ctx.u));
        for k in i..j {
            let ei = cand[k].1 as usize;
            let cs = dot(u2, tab.ax[l.hn][ei].w);
            if !tab.cosbin[(cs.abs() * 4096.0) as usize] { continue; }
            entry(ctx, &w2, o2, u2, l.hn, ei, cs, pieces);
        }
        ctx.desc.pop();
        pieces.pop();
        i = j;
    }
    ctx.bins = bins;
    ctx.cand = cand;
    ctx.scos = scos;
}

fn dfs(ctx: &mut Ctx, w: &Xf, n: usize, depth: usize, any_left: usize, pieces: &mut Vec<(usize, Xf)>) {
    test(ctx, w, n, pieces);
    if depth >= ctx.maxdepth { return; }
    if ctx.join && depth + 1 == ctx.maxdepth { join(ctx, w, n, any_left, pieces); return; }
    let tab = ctx.tab;
    let list = if any_left == 0 { &tab.thru[n] } else { &tab.links[n] };
    for l in list {
        let xl = w.mul(&l.xp);
        if prune(ctx, l.solid, xl.t, pieces) { ctx.pruned += 1; continue; }
        pieces.push((l.solid, xl));
        ctx.desc.push(format!("{}[f{} k{} > f{}]", ctx.lib[l.solid].name, l.g, l.k, l.h));
        dfs(ctx, &w.mul(&l.exit), l.hn, depth + 1, if l.through { any_left } else { any_left - 1 }, pieces);
        ctx.desc.pop();
        pieces.pop();
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut solids_path = "solids.txt".to_string();
    let (mut depth, mut any) = (2usize, 1usize);
    let mut threads = std::thread::available_parallelism().map_or(4, |n| n.get());
    let mut resume = false;
    let (mut maxrad, mut maxpieces) = (f64::MAX, usize::MAX);
    let mut maxfaces = usize::MAX;
    let (mut minorder, mut near) = (24usize, 0usize);
    let mut range: Option<(usize, usize)> = None;
    let mut outdir = "out".to_string();
    let mut only: Option<Vec<String>> = None;
    let mut seeds_only: Option<Vec<String>> = None;
    let mut links_only: Option<Vec<String>> = None;
    let mut use_join = true;
    let mut use_share = true;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--solids" => { solids_path = args[i + 1].clone(); i += 1; }
            "--depth" => { depth = args[i + 1].parse().unwrap(); i += 1; }
            "--any" => { any = args[i + 1].parse().unwrap(); i += 1; }
            "--threads" => { threads = args[i + 1].parse().unwrap(); i += 1; }
            "--min-group" => { minorder = args[i + 1].parse().unwrap(); i += 1; }
            "--near" => { near = args[i + 1].parse().unwrap(); i += 1; }
            "--seed-range" => { let v: Vec<usize> = args[i + 1].split(':').map(|x| x.parse().unwrap()).collect(); range = Some((v[0], v[1])); i += 1; }
            "--resume" => { resume = true; }
            "--no-join" => { use_join = false; }
            "--no-share" => { use_share = false; }
            "--max-radius" => { maxrad = args[i + 1].parse().unwrap(); i += 1; }
            "--max-faces" => { maxfaces = args[i + 1].parse().unwrap(); i += 1; }
            "--max-pieces" => { maxpieces = args[i + 1].parse().unwrap(); i += 1; }
            "--out" => { outdir = args[i + 1].clone(); i += 1; }
            "--links" => { links_only = Some(args[i + 1].split(',').map(|s| s.to_string()).collect()); i += 1; }
            "--only" => { only = Some(args[i + 1].split(',').map(|s| s.to_string()).collect()); i += 1; }
            "--seeds" => { seeds_only = Some(args[i + 1].split(',').map(|s| s.to_string()).collect()); i += 1; }
            x => panic!("unknown arg {}", x),
        }
        i += 1;
    }
    let lib = solids::load(&solids_path);
    let active: Vec<bool> = lib.iter().map(|s| only.as_ref().map_or(true, |o| o.contains(&s.name)) && s.faces.len() <= maxfaces).collect();
    // --links: which solids may be links. "struts" makes the free links prisms and antiprisms
    // (cube and octahedron included) passed base to base, i.e. the straight struts of a
    // node-and-strut frame. With --any K, up to K other pieces may sit in the chain as bends.
    let struts = links_only.as_ref().map_or(false, |l| l.len() == 1 && l[0] == "struts");
    let linkok: Vec<bool> = lib.iter().map(|s| match &links_only {
        None => true,
        Some(_) if struts => true,
        Some(l) => l.contains(&s.name),
    }).collect();
    let tab = build_tables(&lib, &active, &linkok, struts);
    // seeds: (solid, axis, face rep under the axis stabiliser)
    let mut seeds: Vec<(usize, usize, usize)> = Vec::new();
    for (si, s) in lib.iter().enumerate() {
        if !active[si] { continue; }
        if let Some(so) = &seeds_only { if !so.contains(&s.name) { continue; } }
        for (ai, &(u, ord)) in s.axes.iter().enumerate() {
            if !divisors(ord).iter().any(|d| *d >= 3) { continue; }
            let stab: Vec<usize> = (0..s.syms.len()).filter(|&k| dot(mv(&s.syms[k], u), u).abs() > 1.0 - 1e-9).collect();
            for f in 0..s.faces.len() {
                if stab.iter().all(|&k| s.fperm[k][f] >= f) { seeds.push((si, ai, f)); }
            }
        }
    }
    if let Some((a, b)) = range { seeds = seeds[a.min(seeds.len())..b.min(seeds.len())].to_vec(); }
    // Output directory and checkpoint. checkpoint.txt holds the run's settings, every result's
    // dedupe signature and every completed seed, so --resume can skip finished seeds and append.
    std::fs::create_dir_all(format!("{}/off", outdir)).unwrap();
    let config = format!("config depth={} any={} min-group={} near={} only={} seeds={} range={:?} nseeds={}", depth, any, minorder, near,
        only.as_ref().map_or("all".to_string(), |v| v.join(",")), seeds_only.as_ref().map_or("all".to_string(), |v| v.join(",")), range, seeds.len())
        + &(if maxrad < f64::MAX { format!(" max-radius={}", maxrad) } else { String::new() })
        + &links_only.as_ref().map_or(String::new(), |l| format!(" links={}", l.join(",")))
        + &(if maxfaces < usize::MAX { format!(" max-faces={}", maxfaces) } else { String::new() })
        + &(if maxpieces < usize::MAX { format!(" max-pieces={}", maxpieces) } else { String::new() });
    let ckpath = format!("{}/checkpoint.txt", outdir);
    let mut done = vec![false; seeds.len()];
    let mut seen: HashSet<String> = HashSet::new();
    let mut nres = 0usize;
    let resuming = resume && std::path::Path::new(&ckpath).exists();
    if resuming {
        let txt = std::fs::read_to_string(&ckpath).unwrap();
        let mut lines = txt.lines();
        let old = lines.next().unwrap_or("");
        if old != config {
            eprintln!("--resume: {} was made with different settings\n  there: {}\n  now:   {}", outdir, old, config);
            std::process::exit(2);
        }
        for l in lines {
            if let Some(x) = l.strip_prefix("seed ") { if let Ok(i) = x.parse::<usize>() { if i < done.len() { done[i] = true; } } }
            else if let Some(x) = l.strip_prefix("sig ") {
                if let Some((id, sig)) = x.split_once('\t') { if let Ok(id) = id.parse::<usize>() { nres = nres.max(id + 1); seen.insert(sig.to_string()); } }
            }
        }
    }
    let open = |name: &str| std::fs::OpenOptions::new().create(true).write(true).append(resuming).truncate(!resuming).open(format!("{}/{}", outdir, name)).unwrap();
    let mut tsv = open("results.tsv");
    let jsonl = open("results.jsonl");
    let mut ckpt = open("checkpoint.txt");
    *LOG.lock().unwrap() = Some(open("run.log"));
    if !resuming {
        writeln!(tsv, "id\tclass\tgroup\tgenus\tpieces\tA_viol\tQ_missing\tQ_edges\thull\trin\trcirc\troundness\thull_faces\tV\tE\tF\tcontent\tchain").unwrap();
        writeln!(ckpt, "{}", config).unwrap();
    }
    let todo: Vec<usize> = (0..seeds.len()).filter(|i| !done[*i]).collect();
    let nfin = AtomicUsize::new(seeds.len() - todo.len());
    let shared = Mutex::new(Shared { seen, tsv, jsonl, ckpt, n: nres, reasons: Default::default() });
    say(format!("{}{}", if resuming { "resuming: " } else { "" }, config));
    say(format!("through links {:?}", tab.thru.iter().map(|v| v.len()).collect::<Vec<_>>()));
    if use_join { say(format!("join: directions {:?} rotation classes {:?} (through {:?})", tab.ud.iter().map(|v| v.n).collect::<Vec<_>>(),
        tab.rc[1].iter().map(|v| v.len()).collect::<Vec<_>>(), tab.rc[0].iter().map(|v| v.len()).collect::<Vec<_>>())); }
    say(format!("solids {} seeds {} ({} already done) threads {} axtab {:?} links {:?}", lib.len(), seeds.len(), seeds.len() - todo.len(), threads,
        tab.ax.iter().map(|v| v.len()).collect::<Vec<_>>(), tab.links.iter().map(|v| v.len()).collect::<Vec<_>>()));
    let next = AtomicUsize::new(0);
    let (tests, hits) = (AtomicU64::new(0), AtomicU64::new(0));
    let t0 = std::time::Instant::now();
    std::thread::scope(|sc| {
        if std::env::var("WATCH").is_ok() {
            sc.spawn(|| while G_DONE.load(Ordering::Relaxed) == 0 {
                // sleep in short steps so the process exits promptly when the work is done
                for _ in 0..40 { if G_DONE.load(Ordering::Relaxed) != 0 { return; } std::thread::sleep(std::time::Duration::from_millis(250)); }
                say(format!("watch {:.0}s seeds {}/{} tests {} hits {} quick {} validate_cpu_s {:.1}", t0.elapsed().as_secs_f64(), nfin.load(Ordering::Relaxed), seeds.len(), G_TESTS.load(Ordering::Relaxed),
                    G_HITS.load(Ordering::Relaxed), G_QUICK.load(Ordering::Relaxed), G_VAL_NS.load(Ordering::Relaxed) as f64 / 1e9));
                let parts: Vec<String> = (0..12).filter(|k| G_R[*k].load(Ordering::Relaxed) > 0).map(|k| format!("{} n={} cpu={:.1}s", RNAMES[k], G_R[k].load(Ordering::Relaxed), G_RNS[k].load(Ordering::Relaxed) as f64 / 1e9)).collect();
                say(format!("  {}", parts.join(" | ")));
            });
        }
        let workers: Vec<_> = (0..threads).map(|_| {
            sc.spawn(|| loop {
                let k = next.fetch_add(1, Ordering::Relaxed);
                if k >= todo.len() { break; }
                let i = todo[k];
                let (si, ai, f) = seeds[i];
                let s = &lib[si];
                let u = s.axes[ai].0;
                let mut inside = (f64::MIN, f64::MAX);
                for f in 0..s.faces.len() {
                    let (a, b) = (dot(s.fnorm[f], u), dot(s.fnorm[f], s.v[s.faces[f][0]]) - 1e-6);
                    if a > 0.0 { inside.1 = inside.1.min(b / a); } else if a < 0.0 { inside.0 = inside.0.max(b / a); }
                }
                let mut ctx = Ctx { lib: &lib, tab: &tab, a: si, u: s.axes[ai].0, maxdepth: depth, ordered: seeds_only.is_none(), minorder, near, shared: &shared, outdir: &outdir, tests: 0, hits: 0,
                    desc: vec![format!("{}[ax{} f{}]", s.name, s.axes[ai].1, f)], reasons: Default::default(), scratch: Default::default(),
                    pdiv: divisors(s.axes[ai].1).into_iter().filter(|x| *x >= 3).collect(), inside, csbuf: Vec::new(), st: Default::default(),
                    axrot: Vec::new(), maxrad, pruned: 0, join: use_join, sbin: Box::new([false; 65536]), share: use_share, check: std::env::var("SHARE_CHECK").is_ok(), gc: Default::default(), gce: Vec::new(), gcg: Vec::new(), gc_np: 1, ex: [0.0; 3], ey: [0.0; 3], pre_ov: 0, pre_q: 0, bins: Vec::new(), scos: Vec::new(), cand: Vec::new() };
                for x in &tab.pq {
                    if !ctx.pdiv.contains(&x.0) || x.3.len() < minorder { continue; }
                    let b = ((x.2 + 1.0) * 4096.0) as usize;
                    for k in b.saturating_sub(1)..=b + 1 { ctx.sbin[k] = true; }
                    if !ctx.scos.iter().any(|c| (c - x.2).abs() < 1e-9) { ctx.scos.push(x.2); }
                }
                let a = if u[0].abs() < 0.6 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
                ctx.ex = unit(cross(u, a));
                ctx.ey = cross(u, ctx.ex);
                ctx.axrot = ctx.pdiv.iter().map(|&p| rot_sc(u, tab.sc[p])).collect();
                ctx.scratch.max_pieces = maxpieces;
                let mut pieces = vec![(ctx.a, Xf::ident())];
                let w = s.frames[f][0];
                dfs(&mut ctx, &w, s.faces[f].len(), 0, any, &mut pieces);
                ctx.st.flush();
                if ctx.pre_ov > 0 { *ctx.reasons.entry("prefix-overlap").or_insert(0) += ctx.pre_ov; }
                if ctx.pre_q > 0 { *ctx.reasons.entry("prefix-q").or_insert(0) += ctx.pre_q; }
                if ctx.pruned > 0 { *ctx.reasons.entry("pruned-subtrees").or_insert(0) += ctx.pruned; }
                if ctx.st.quick_total > 0 { *ctx.reasons.entry("overlap-quick").or_insert(0) += ctx.st.quick_total; }
                tests.fetch_add(ctx.tests, Ordering::Relaxed);
                {
                    let mut sh = shared.lock().unwrap();
                    for (k, v) in &ctx.reasons { *sh.reasons.entry(*k).or_insert(0) += *v; }
                    // the seed's results are already on disk, so it is safe to mark it complete
                    writeln!(sh.ckpt, "seed {}", i).unwrap();
                }
                hits.fetch_add(ctx.hits, Ordering::Relaxed);
                let nf = nfin.fetch_add(1, Ordering::Relaxed) + 1;
                if nf % 10 == 0 || nf == seeds.len() { say(format!("seeds done {}/{} {:.0}s", nf, seeds.len(), t0.elapsed().as_secs_f64())); }
            })
        }).collect();
        for w in workers { w.join().unwrap(); }
        G_DONE.store(1, Ordering::Relaxed);
    });
    let sh = shared.lock().unwrap();
    say(format!("done {:.1}s tests {} hits {} unique {} reasons {:?}{}", t0.elapsed().as_secs_f64(), tests.load(Ordering::Relaxed), hits.load(Ordering::Relaxed), sh.n, sh.reasons,
        if todo.len() < seeds.len() { " (counts cover this session only; unique is the total)" } else { "" }));
}
