//! Grow a known toroid outward: attach further solids face to face to its pieces, take their
//! orbit under the toroid's own group, and keep what is still an (R)(A)(Q)(T) toroid on a
//! convex hull other than the one it started from.
//!
//! Unlike a search hit, a decoration does not have to end on a symmetry axis, so this reaches
//! shapes the chain search cannot: a ring of pieces around a node, a piece in a window, a cap.
//! Level 1 tries every solid on every face of every representative piece. Each further level
//! adds either another level-1 decoration or a solid on a face of a piece added before. The
//! last level runs with the (Q) screens on; earlier levels need only be free of collisions.
use crate::geom::*;
use crate::solids::Solid;
use crate::validate;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

const FLIP: M3 = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]];

fn nums(s: &str) -> Vec<f64> { s.split(|c: char| c == ',' || c == '[' || c == ']' || c == '}').filter(|x| !x.trim().is_empty()).map(|x| x.trim().parse().unwrap()).collect() }
fn xf_of(v: &[f64]) -> Xf { Xf { r: [[v[0], v[1], v[2]], [v[3], v[4], v[5]], [v[6], v[7], v[8]]], t: [v[9], v[10], v[11]] } }

/// Every placement of every solid on every face of `parent`.
fn attach(lib: &[Solid], active: &[bool], parent: &(usize, Xf), out: &mut Vec<(usize, Xf)>) {
    let ps = &lib[parent.0];
    for h in 0..ps.faces.len() {
        let n = ps.faces[h].len();
        let w = parent.1.mul(&ps.frames[h][0]);
        for (t, s) in lib.iter().enumerate() {
            if !active[t] { continue; }
            for &(g, k) in &s.flags {
                if s.faces[g].len() != n { continue; }
                let fi = s.frames[g][k].inv();
                out.push((t, w.mul(&Xf { r: mmul(&FLIP, &fi.r), t: mv(&FLIP, fi.t) })));
            }
        }
    }
}

struct Out { hulls: HashSet<String>, seen: HashSet<String>, tsv: std::fs::File, jsonl: std::fs::File, n: usize, reasons: HashMap<&'static str, u64> }

struct Cx<'a> { outward: bool, lib: &'a [Solid], active: &'a [bool], grp: &'a [Xf], centre: V3, gens: Vec<Xf>, nbase: usize, rin0: f64, maxrad: f64, depth: usize, out: &'a str, shared: &'a Mutex<Out>, first: &'a [(usize, Xf)] }

fn sig(r: &validate::Res, centre: V3) -> String {
    let (mut m1, mut m3) = (0.0, 0.0);
    for p in &r.pts { let d = dist(*p, centre); m1 += d; m3 += d * d * d; }
    format!("{} {} {} {} {} {:.4} {:.3} {:.2} {:?}", r.gorder, r.genus, r.nv, r.ne, r.nf, r.rcirc, m1, m3, r.counts)
}

/// The hull as `--new-hulls` files spell it, but with each face size listed once: how often a
/// size is listed depends on which pieces are the representatives, not only on the hull.
fn hull_key(r: &validate::Res) -> String {
    let hf: Vec<String> = r.hull_faces.iter().map(|n| n.to_string()).collect();
    reduce(&format!("{}\t{:.4}\t{:.4}\t{}", r.gorder, r.rcirc, r.rin, hf.join(",")))
}

fn reduce(key: &str) -> String {
    let f: Vec<&str> = key.split('\t').collect();
    if f.len() < 4 { return key.to_string(); }
    let mut sizes: Vec<usize> = f[3].split(',').filter_map(|x| x.parse().ok()).collect();
    sizes.sort();
    sizes.dedup();
    format!("{}\t{}\t{}\t{}", f[0], f[1], f[2], sizes.iter().map(|n| n.to_string()).collect::<Vec<_>>().join(","))
}

/// Validate `reps` (the base plus decorations). Returns the result when the union is free of
/// collisions and a manifold, i.e. something a further level may build on.
fn try_one(cx: &Cx, sc: &mut validate::Scratch, reasons: &mut HashMap<&'static str, u64>, reps: &[(usize, Xf)], last: bool) -> Option<validate::Res> {
    let lib = cx.lib;
    let (s, x) = reps[reps.len() - 1];
    if dist(x.t, cx.centre) + lib[s].radius > cx.maxrad { *reasons.entry("too-far").or_insert(0) += 1; return None; }
    // the same solid already sits there (validate would quietly merge the two)
    for (s2, x2) in &reps[..reps.len() - 1] {
        if *s2 == s && cx.grp.iter().any(|g| dist(g.ap(x2.t), x.t) < 1e-5) { *reasons.entry("duplicate").or_insert(0) += 1; return None; }
    }
    // --outward: every added piece must reach the old hull
    if cx.outward && dist(x.t, cx.centre) + lib[s].radius < cx.rin0 - 1e-6 { *reasons.entry("inside").or_insert(0) += 1; return None; }
    // nothing added reaches the old hull: the hull cannot have changed
    if last && reps[cx.nbase..].iter().all(|(s, x)| dist(x.t, cx.centre) + lib[*s].radius < cx.rin0 - 1e-6) { *reasons.entry("inside").or_insert(0) += 1; return None; }
    let r = match validate::validate(lib, reps, cx.centre, last, sc) {
        Err(e) => { *reasons.entry(e).or_insert(0) += 1; return None; }
        Ok(r) => r,
    };
    *reasons.entry("valid").or_insert(0) += 1;
    if r.a_viol == 0 && r.q_ok {
        *reasons.entry("RAQT").or_insert(0) += 1;
        let hull = hull_key(&r);
        let mut sh = cx.shared.lock().unwrap();
        if !sh.hulls.insert(hull) { *reasons.entry("known-hull").or_insert(0) += 1; }
        else if sh.seen.insert(sig(&r, cx.centre)) {
            let id = sh.n;
            sh.n += 1;
            let names: Vec<String> = r.counts.iter().map(|(s, c)| format!("{}x{}", c, lib[*s].name)).collect();
            let hf: Vec<String> = r.hull_faces.iter().map(|n| n.to_string()).collect();
            let chain: Vec<&str> = reps[cx.nbase..].iter().map(|(s, _)| lib[*s].name.as_str()).collect();
            writeln!(sh.tsv, "{}\tRAQT\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\tbase + {}",
                id, r.gorder, r.genus, r.npieces, r.a_viol, r.q_missing, r.q_edges, if r.hull_regular { "regular" } else { "irregular" }, r.rin, r.rcirc,
                if r.rin < 1e9 { r.rin / r.rcirc } else { 0.0 }, hf.join(","), r.nv, r.ne, r.nf, names.join("+"), chain.join(" + ")).unwrap();
            let rj: Vec<String> = reps.iter().map(|(s, x)| format!("{{\"solid\":\"{}\",\"xf\":{}}}", lib[*s].name, crate::xf_json(x))).collect();
            let gj: Vec<String> = cx.gens.iter().map(crate::xf_json).collect();
            writeln!(sh.jsonl, "{{\"id\":{},\"good\":true,\"group\":{},\"genus\":{},\"centre\":[{:.15e},{:.15e},{:.15e}],\"reps\":[{}],\"gens\":[{}]}}",
                id, r.gorder, r.genus, cx.centre[0], cx.centre[1], cx.centre[2], rj.join(","), gj.join(",")).unwrap();
            crate::write_off(&format!("{}/off/{:06}.off", cx.out, id), &r, lib);
            eprintln!("new hull: R {:.4} genus {} pieces {} hull faces {} : base + {}", r.rcirc, r.genus, r.npieces, hf.join(","), chain.join(" + "));
        }
    }
    Some(r)
}

/// `reps` holds the base and the decorations so far; `from` is the first level-1 decoration
/// that may still be added (they are added in list order so that no pair is tried twice).
fn grow(cx: &Cx, sc: &mut validate::Scratch, reasons: &mut HashMap<&'static str, u64>, reps: &mut Vec<(usize, Xf)>, from: usize, level: usize) {
    let last = level == cx.depth;
    let mut cands: Vec<(usize, Xf, usize)> = Vec::new();
    let mut tmp = Vec::new();
    for p in &reps[cx.nbase..] { attach(cx.lib, cx.active, p, &mut tmp); }
    for (s, x) in tmp { cands.push((s, x, from)); }
    for i in from..cx.first.len() { cands.push((cx.first[i].0, cx.first[i].1, i + 1)); }
    for (s, x, nf) in cands {
        reps.push((s, x));
        if try_one(cx, sc, reasons, reps, last).is_some() && !last { grow(cx, sc, reasons, reps, nf, level + 1); }
        reps.pop();
    }
}

pub fn decorate(lib: &[Solid], active: &[bool], file: &str, id: usize, depth: usize, maxpieces: usize, threads: usize, hulls: Option<HashSet<String>>, out: &str, maxrad: f64, beam: usize, outward: bool) {
    let text = std::fs::read_to_string(file).expect("cannot read the results.jsonl to decorate");
    let line = text.lines().find(|l| l.starts_with(&format!("{{\"id\":{},", id))).expect("no result with that id");
    let after = |key: &str| -> &str { &line[line.find(key).unwrap() + key.len()..] };
    let c = { let s = after("\"centre\":["); nums(&s[..s.find(']').unwrap()]) };
    let centre = [c[0], c[1], c[2]];
    let mut base: Vec<(usize, Xf)> = Vec::new();
    for part in line.split("{\"solid\":\"").skip(1) {
        let name = &part[..part.find('"').unwrap()];
        let a = part.find("\"xf\":[").unwrap() + 6;
        let b = a + part[a..].find(']').unwrap();
        base.push((lib.iter().position(|s| s.name == name).expect("unknown solid"), xf_of(&nums(&part[a..b]))));
    }
    let g = nums(after("\"gens\":["));
    let gens = vec![xf_of(&g[0..12]), xf_of(&g[12..24])];
    let grp = closure(&gens, 120).unwrap();
    let new_scratch = || { let mut sc = validate::Scratch::default(); sc.max_pieces = maxpieces; sc.grp.extend_from_slice(&grp); sc };
    let mut sc = new_scratch();
    let r0 = validate::validate(lib, &base, centre, false, &mut sc).expect("the base does not validate");
    let hull0 = hull_key(&r0);
    eprintln!("base: {} pieces genus {} R {:.4} rin {:.4} group {}", r0.npieces, r0.genus, r0.rcirc, r0.rin, r0.gorder);
    std::fs::create_dir_all(format!("{}/off", out)).unwrap();
    let mut tsv = std::fs::File::create(format!("{}/results.tsv", out)).unwrap();
    writeln!(tsv, "id\tclass\tgroup\tgenus\tpieces\tA_viol\tQ_missing\tQ_edges\thull\trin\trcirc\troundness\thull_faces\tV\tE\tF\tcontent\tchain").unwrap();
    let jsonl = std::fs::File::create(format!("{}/results.jsonl", out)).unwrap();
    let mut hs: HashSet<String> = hulls.unwrap_or_default().iter().map(|h| reduce(h)).collect();
    hs.insert(hull0);
    let shared = Mutex::new(Out { hulls: hs, seen: HashSet::new(), tsv, jsonl, n: 0, reasons: HashMap::new() });
    let t0 = std::time::Instant::now();
    // level 1: every solid on every face of every base piece
    let mut c1: Vec<(usize, Xf)> = Vec::new();
    for p in &base { attach(lib, active, p, &mut c1); }
    let first: Mutex<Vec<(usize, Xf, String)>> = Mutex::new(Vec::new());
    let next = AtomicUsize::new(0);
    let run = |work: &(dyn Fn(usize, &mut validate::Scratch, &mut HashMap<&'static str, u64>) + Sync), n: usize| {
        next.store(0, Ordering::Relaxed);
        std::thread::scope(|scope| {
            for _ in 0..threads {
                scope.spawn(|| {
                    let mut sc = new_scratch();
                    let mut reasons: HashMap<&'static str, u64> = HashMap::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= n { break; }
                        work(i, &mut sc, &mut reasons);
                    }
                    let mut sh = shared.lock().unwrap();
                    for (k, v) in reasons { *sh.reasons.entry(k).or_insert(0) += v; }
                });
            }
        });
    };
    {
        let cx = Cx { outward, lib, active, grp: &grp, centre, gens: gens.clone(), nbase: base.len(), rin0: r0.rin, maxrad, depth, out, shared: &shared, first: &[] };
        run(&|i, sc, reasons| {
            let mut reps = base.clone();
            reps.push(c1[i]);
            // collision test without the (Q) screens, so that what survives can be built on
            if let Some(r) = try_one(&cx, sc, reasons, &reps, false) { first.lock().unwrap().push((c1[i].0, c1[i].1, sig(&r, centre))); }
        }, c1.len());
    }
    let mut fv = first.into_inner().unwrap();
    fv.sort_by(|a, b| a.2.cmp(&b.2).then(a.0.cmp(&b.0)));
    fv.dedup_by(|a, b| a.2 == b.2);
    let first: Vec<(usize, Xf)> = fv.into_iter().map(|(s, x, _)| (s, x)).collect();
    eprintln!("level 1: {} placements tried, {} distinct collision-free decorations, {:.1}s", c1.len(), first.len(), t0.elapsed().as_secs_f64());
    if beam > 0 {
        // Beam search: every state is the base plus some decorations that reach beyond the old
        // hull. A state is scored by the share of the hull edges found through its representative
        // pieces that are not edges of the toroid (0 would be quasi-convex, 1 if none was found);
        // the `beam` best distinct states of a level are each extended by one more decoration,
        // on a base piece or on a piece added before. A heuristic: the score is a lower bound.
        let cx = Cx { outward, lib, active, grp: &grp, centre, gens: gens.clone(), nbase: base.len(), rin0: r0.rin, maxrad, depth, out, shared: &shared, first: &first };
        let outward = |c: &(usize, Xf)| dist(c.1.t, centre) + lib[c.0].radius > r0.rin + 1e-6;
        let c1: Vec<(usize, Xf)> = c1.iter().filter(|c| outward(c)).cloned().collect();
        let mut states: Vec<Vec<(usize, Xf)>> = vec![base.clone()];
        let mut seen: HashSet<String> = HashSet::new();
        for level in 1..=depth {
            let found: Mutex<Vec<(f64, usize, String, Vec<(usize, Xf)>)>> = Mutex::new(Vec::new());
            run(&|i, sc, reasons| {
                let st = &states[i];
                let mut cands = c1.clone();
                for p in &st[base.len()..] { attach(lib, active, p, &mut cands); }
                let mut reps = st.clone();
                let mut loc = Vec::new();
                for c in cands {
                    if !outward(&c) { continue; }
                    reps.push(c);
                    if let Some(r) = try_one(&cx, sc, reasons, &reps, false) {
                        if r.rcirc > r0.rcirc + 1e-6 { loc.push((if r.q_edges == 0 { 1.0 } else { r.q_missing as f64 / r.q_edges as f64 }, r.npieces, sig(&r, centre), reps.clone())); }
                    }
                    reps.pop();
                }
                found.lock().unwrap().append(&mut loc);
            }, states.len());
            let mut fv = found.into_inner().unwrap();
            let total = fv.len();
            fv.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
            fv.dedup_by(|a, b| a.2 == b.2);
            fv.retain(|f| seen.insert(f.2.clone()));
            fv.truncate(beam);
            if fv.is_empty() { eprintln!("level {}: nothing to extend", level); break; }
            let names = |v: &Vec<(usize, Xf)>| v[base.len()..].iter().map(|(s, _)| lib[*s].name.as_str()).collect::<Vec<_>>().join(" + ");
            eprintln!("level {}: {} states from {}, {} kept, hull-edge defect best {:.3} worst kept {:.3}, {:.0}s; best: {}", level, total, states.len(), fv.len(), fv[0].0, fv[fv.len() - 1].0, t0.elapsed().as_secs_f64(), names(&fv[0].3));
            states = fv.into_iter().map(|f| f.3).collect();
        }
    } else if depth >= 2 {
        let cx = Cx { outward, lib, active, grp: &grp, centre, gens: gens.clone(), nbase: base.len(), rin0: r0.rin, maxrad, depth, out, shared: &shared, first: &first };
        let done = AtomicUsize::new(0);
        run(&|i, sc, reasons| {
            let mut reps = base.clone();
            reps.push(first[i]);
            grow(&cx, sc, reasons, &mut reps, i + 1, 2);
            let d = done.fetch_add(1, Ordering::Relaxed) + 1;
            if d % 200 == 0 { eprintln!("level-1 decorations expanded {}/{} {:.0}s", d, first.len(), t0.elapsed().as_secs_f64()); }
        }, first.len());
    }
    let sh = shared.lock().unwrap();
    let mut rs: Vec<_> = sh.reasons.iter().collect();
    rs.sort();
    println!("done {:.1}s new hulls {} reasons {:?}", t0.elapsed().as_secs_f64(), sh.n, rs);
}
