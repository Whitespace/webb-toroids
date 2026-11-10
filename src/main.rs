mod geom;
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

struct AxEntry { solid: usize, face: usize, k: usize, order: usize, xp: Xf, p: V3, w: V3 }
#[derive(Clone)]
struct Link { solid: usize, g: usize, k: usize, h: usize, hn: usize, xp: Xf, exit: Xf, through: bool }

struct Tables { ax: Vec<Vec<AxEntry>>, links: Vec<Vec<Link>>, thru: Vec<Vec<Link>>, cosset: Vec<f64>, pq: Vec<(usize, usize, f64, Vec<(usize, usize)>)>, cosbin: Vec<bool> }

fn build_tables(lib: &[Solid], active: &[bool]) -> Tables {
    let mut ax: Vec<Vec<AxEntry>> = (0..11).map(|_| Vec::new()).collect();
    let mut links: Vec<Vec<Link>> = (0..11).map(|_| Vec::new()).collect();
    for (si, s) in lib.iter().enumerate() {
        if !active[si] { continue; }
        for &(g, k) in &s.flags {
            let n = s.faces[g].len();
            let fi = s.frames[g][k].inv();
            let xp = Xf { r: mmul(&FLIP, &fi.r), t: mv(&FLIP, fi.t) };
            for &(u, order) in &s.axes {
                ax[n].push(AxEntry { solid: si, face: g, k, order, xp, p: xp.t, w: mv(&xp.r, u) });
            }
            for h in 0..s.faces.len() {
                if h == g { continue; }
                let through = dot(s.fnorm[h], s.fnorm[g]) < -1.0 + 1e-9;
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
    Tables { ax, links, thru, cosset, pq, cosbin }
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
    na: usize,
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
    scratch: Vec<V3>,
}

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
    let wi = w.inv();
    let o = wi.t;
    let u = mv(&wi.r, ctx.u);
    let tab = ctx.tab;
    ctx.tests += tab.ax[n].len() as u64;
    G_TESTS.fetch_add(tab.ax[n].len() as u64, Ordering::Relaxed);
    for e in &tab.ax[n] {
        let cs = dot(u, e.w);
        let c = cs.abs();
        if !tab.cosbin[(c * 4096.0) as usize] { continue; }
        if !tab.cosset.iter().any(|x| (x - c).abs() < 1e-7) { continue; }
        let nn = cross(u, e.w);
        let d = sub(e.p, o);
        if dot(d, nn).abs() > 1e-7 { continue; }
        let bq: Vec<usize> = { let d: Vec<usize> = divisors(e.order).into_iter().filter(|x| *x >= 3).collect(); if d.is_empty() { vec![2] } else { d } };
        if ctx.ordered && bq[0] >= 3 && e.solid < ctx.a { continue; }
        let s = dot(cross(d, e.w), nn) / dot(nn, nn);
        let centre = scale(ctx.u, s);
        let xb = w.mul(&e.xp);
        // a piece sitting on the centre gives a star of arms, not a shell; skip those
        if s.abs() < 1e-6 || dist(centre, xb.t) < 1e-6 { continue; }
        // centre strictly inside the seed piece: it would meet its own images
        let sa = &ctx.lib[ctx.a];
        if (0..sa.faces.len()).all(|f| dot(sa.fnorm[f], sub(centre, sa.v[sa.faces[f][0]])) < -1e-6) { continue; }
        let wb = mv(&w.r, e.w);
        pieces.push((e.solid, xb));
        for p in divisors(ctx.na).into_iter().filter(|x| *x >= 3) {
            for &q in &bq {
                let rec = match tab.pq.iter().find(|x| x.0 == p && x.1 == q && (x.2 - cs).abs() < 1e-7 && x.3.len() >= ctx.minorder) { Some(x) => &x.3, None => continue };
                ctx.hits += 1;
                G_HITS.fetch_add(1, Ordering::Relaxed);
                let gens = [Xf::about(centre, ctx.u, 2.0 * std::f64::consts::PI / p as f64), Xf::about(centre, wb, 2.0 * std::f64::consts::PI / q as f64)];
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
                if quick { G_QUICK.fetch_add(1, Ordering::Relaxed); *ctx.reasons.entry("overlap-quick").or_insert(0) += 1; continue; }
                let t0 = std::time::Instant::now();
                let grp = from_recipe(&gens, rec);
                let res = validate::validate(ctx.lib, pieces, &grp, centre, ctx.near == 0, &mut ctx.scratch);
                let ns = t0.elapsed().as_nanos() as u64;
                G_VAL_NS.fetch_add(ns, Ordering::Relaxed);
                if let Err(r) = &res { let k = RNAMES.iter().position(|x| x == r).unwrap_or(11); G_R[k].fetch_add(1, Ordering::Relaxed); G_RNS[k].fetch_add(ns, Ordering::Relaxed); }
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

fn dfs(ctx: &mut Ctx, w: &Xf, n: usize, depth: usize, any_left: usize, pieces: &mut Vec<(usize, Xf)>) {
    test(ctx, w, n, pieces);
    if depth >= ctx.maxdepth { return; }
    let tab = ctx.tab;
    let list = if any_left == 0 { &tab.thru[n] } else { &tab.links[n] };
    for l in list {
        pieces.push((l.solid, w.mul(&l.xp)));
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
    let (mut minorder, mut near) = (24usize, 0usize);
    let mut range: Option<(usize, usize)> = None;
    let mut outdir = "out".to_string();
    let mut only: Option<Vec<String>> = None;
    let mut seeds_only: Option<Vec<String>> = None;
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
            "--out" => { outdir = args[i + 1].clone(); i += 1; }
            "--only" => { only = Some(args[i + 1].split(',').map(|s| s.to_string()).collect()); i += 1; }
            "--seeds" => { seeds_only = Some(args[i + 1].split(',').map(|s| s.to_string()).collect()); i += 1; }
            x => panic!("unknown arg {}", x),
        }
        i += 1;
    }
    let lib = solids::load(&solids_path);
    let active: Vec<bool> = lib.iter().map(|s| only.as_ref().map_or(true, |o| o.contains(&s.name))).collect();
    let tab = build_tables(&lib, &active);
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
        only.as_ref().map_or("all".to_string(), |v| v.join(",")), seeds_only.as_ref().map_or("all".to_string(), |v| v.join(",")), range, seeds.len());
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
                let mut ctx = Ctx { lib: &lib, tab: &tab, a: si, u: s.axes[ai].0, na: s.axes[ai].1, maxdepth: depth, ordered: seeds_only.is_none(), minorder, near, shared: &shared, outdir: &outdir, tests: 0, hits: 0,
                    desc: vec![format!("{}[ax{} f{}]", s.name, s.axes[ai].1, f)], reasons: Default::default(), scratch: Vec::new() };
                let mut pieces = vec![(ctx.a, Xf::ident())];
                let w = s.frames[f][0];
                dfs(&mut ctx, &w, s.faces[f].len(), 0, any, &mut pieces);
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
