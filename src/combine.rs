//! Two strut orbits: unions of two collision-free chains that share the seed and the group.
//!
//! A search run with `--partials` writes every hit whose orbit is free of collisions to
//! `partials.tsv`, whether or not it is a toroid by itself, keyed by where the seed sits in
//! the group (the key of the chain-prefix cache: group, centre, azimuth of a reference axis).
//! Two chains with the same seed, seed axis and key live in the same group around the same seed,
//! so their union is a candidate with two chain orbits. Each union goes through the same
//! `validate` as a search hit.
use crate::geom::*;
use crate::solids::Solid;
use crate::validate;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

struct Part { cls: String, centre: V3, gens: [Xf; 2], reps: Vec<(usize, Xf)>, desc: String }

fn xf_of(v: &[f64]) -> Xf { Xf { r: [[v[0], v[1], v[2]], [v[3], v[4], v[5]], [v[6], v[7], v[8]]], t: [v[9], v[10], v[11]] } }

pub fn xf_txt(x: &Xf) -> String {
    let v = [x.r[0][0], x.r[0][1], x.r[0][2], x.r[1][0], x.r[1][1], x.r[1][2], x.r[2][0], x.r[2][1], x.r[2][2], x.t[0], x.t[1], x.t[2]];
    v.iter().map(|a| format!("{:.17e}", a)).collect::<Vec<_>>().join(" ")
}

fn nums(s: &str) -> Vec<f64> { s.split(' ').filter(|x| !x.is_empty()).map(|x| x.parse().unwrap()).collect() }

fn pkey(s: usize, c: V3) -> (usize, [i64; 3]) { (s, [(c[0] * 1e4).round() as i64, (c[1] * 1e4).round() as i64, (c[2] * 1e4).round() as i64]) }

struct Out { hulls: Option<HashSet<String>>, seen: HashSet<String>, tsv: std::fs::File, jsonl: std::fs::File, n: usize, reasons: HashMap<&'static str, u64> }

/// `base_raqt`: only unions in which at least one chain is a toroid by itself (growing a known one).
/// `rounds`, `beam`: how many chain orbits to add in turn, and how many unions survive each round.
pub fn combine(lib: &[Solid], dir: &str, maxpieces: usize, threads: usize, maxpairs: usize, base_raqt: bool, rounds: usize, beam: usize, hulls: Option<HashSet<String>>) {
    let text = std::fs::read_to_string(format!("{}/partials.tsv", dir)).expect("no partials.tsv in that run; search with --partials first");
    let mut buckets: HashMap<(usize, usize, i64, i64, i64), Vec<Part>> = HashMap::new();
    let mut lines: HashSet<&str> = HashSet::new();
    let mut nparts = 0usize;
    for line in text.lines() {
        if !lines.insert(line) { continue; }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 { continue; }
        let key = (f[0].parse().unwrap(), f[1].parse().unwrap(), f[2].parse().unwrap(), f[3].parse().unwrap(), f[4].parse().unwrap());
        let c = nums(f[6]);
        let g = nums(f[7]);
        let reps: Vec<(usize, Xf)> = f[8].split(';').map(|r| { let v = nums(r); (v[0] as usize, xf_of(&v[1..13])) }).collect();
        buckets.entry(key).or_default().push(Part { cls: f[5].to_string(), centre: [c[0], c[1], c[2]], gens: [xf_of(&g[0..12]), xf_of(&g[12..24])], reps, desc: f[9].to_string() });
        nparts += 1;
    }
    drop(lines);
    let mut bl: Vec<Vec<Part>> = buckets.into_values().filter(|b| b.len() > 1).collect();
    bl.sort_by_key(|b| std::cmp::Reverse(b.len()));
    let pairs: usize = bl.iter().map(|b| b.len() * (b.len() - 1) / 2).sum();
    eprintln!("partials {} placements with more than one {} largest {} pairs before dedupe {}", nparts, bl.len(), bl.first().map_or(0, |b| b.len()), pairs);
    let out = format!("{}/combined", dir);
    std::fs::create_dir_all(format!("{}/off", out)).unwrap();
    let mut tsv = std::fs::File::create(format!("{}/results.tsv", out)).unwrap();
    writeln!(tsv, "id\tclass\tgroup\tgenus\tpieces\tA_viol\tQ_missing\tQ_edges\thull\trin\trcirc\troundness\thull_faces\tV\tE\tF\tcontent\tchain").unwrap();
    let jsonl = std::fs::File::create(format!("{}/results.jsonl", out)).unwrap();
    let shared = Mutex::new(Out { hulls, seen: HashSet::new(), tsv, jsonl, n: 0, reasons: HashMap::new() });
    let next = AtomicUsize::new(0);
    let t0 = std::time::Instant::now();
    std::thread::scope(|sc| {
        for _ in 0..threads {
            sc.spawn(|| {
                let mut scratch = validate::Scratch::default();
                scratch.max_pieces = maxpieces;
                let mut reasons: HashMap<&'static str, u64> = HashMap::new();
                loop {
                    let bi = next.fetch_add(1, Ordering::Relaxed);
                    if bi >= bl.len() { break; }
                    let b = &bl[bi];
                    let grp = closure(&b[0].gens, 120).unwrap();
                    // each chain's orbit as a set of (solid, centre); chains with the same orbit are one
                    let mut sets: Vec<HashSet<(usize, [i64; 3])>> = Vec::new();
                    let mut keep: Vec<usize> = Vec::new();
                    let mut sigs: HashSet<Vec<(usize, [i64; 3])>> = HashSet::new();
                    for (i, p) in b.iter().enumerate() {
                        let mut v: Vec<(usize, [i64; 3])> = Vec::new();
                        for (s, x) in &p.reps[1..] { for g in &grp { v.push(pkey(*s, g.ap(x.t))); } }
                        v.sort();
                        v.dedup();
                        if !sigs.insert(v.clone()) { continue; }
                        keep.push(i);
                        sets.push(v.into_iter().collect());
                    }
                    scratch.grp.clear();
                    scratch.grp.extend_from_slice(&grp);
                    // validate a union and, if it is a toroid, record it; returns (genus, hull key)
                    let mut try_union = |reps: &[(usize, Xf)], chain: &str, reasons: &mut HashMap<&'static str, u64>, record: &dyn Fn(i64, &str) -> bool| -> Option<(i64, String)> {
                        match validate::validate(lib, reps, b[0].centre, true, &mut scratch) {
                            Err(e) => { *reasons.entry(e).or_insert(0) += 1; None }
                            Ok(r) => {
                                *reasons.entry("valid").or_insert(0) += 1;
                                if !(r.a_viol == 0 && r.q_ok) { return None; }
                                *reasons.entry("RAQT").or_insert(0) += 1;
                                let hf: Vec<String> = r.hull_faces.iter().map(|n| n.to_string()).collect();
                                let hull = format!("{}\t{:.4}\t{:.4}\t{}", r.gorder, r.rcirc, r.rin, hf.join(","));
                                if !record(r.genus, &hull) { return Some((r.genus, hull)); }
                                let names: Vec<String> = r.counts.iter().map(|(s, c)| format!("{}x{}", c, lib[*s].name)).collect();
                                let (mut m1, mut m3) = (0.0, 0.0);
                                for p in &r.pts { let d = dist(*p, b[0].centre); m1 += d; m3 += d * d * d; }
                                let sig = format!("{} {} {} {} {} {:.4} {:.3} {:.2} {}", r.gorder, r.genus, r.nv, r.ne, r.nf, r.rcirc, m1, m3, names.join("+"));
                                let mut sh = shared.lock().unwrap();
                                if !sh.seen.insert(sig) { return Some((r.genus, hull)); }
                                let id = sh.n;
                                sh.n += 1;
                                writeln!(sh.tsv, "{}\tRAQT\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}",
                                    id, r.gorder, r.genus, r.npieces, r.a_viol, r.q_missing, r.q_edges, if r.hull_regular { "regular" } else { "irregular" }, r.rin, r.rcirc,
                                    if r.rin < 1e9 { r.rin / r.rcirc } else { 0.0 }, hf.join(","), r.nv, r.ne, r.nf, names.join("+"), chain).unwrap();
                                let rj: Vec<String> = reps.iter().map(|(s, x)| format!("{{\"solid\":\"{}\",\"xf\":{}}}", lib[*s].name, crate::xf_json(x))).collect();
                                let gj: Vec<String> = b[0].gens.iter().map(crate::xf_json).collect();
                                writeln!(sh.jsonl, "{{\"id\":{},\"good\":true,\"group\":{},\"genus\":{},\"centre\":[{:.15e},{:.15e},{:.15e}],\"reps\":[{}],\"gens\":[{}]}}",
                                    id, r.gorder, r.genus, b[0].centre[0], b[0].centre[1], b[0].centre[2], rj.join(","), gj.join(",")).unwrap();
                                crate::write_off(&format!("{}/off/{:06}.off", out, id), &r, lib);
                                Some((r.genus, hull))
                            }
                        }
                    };
                    if base_raqt {
                        // Grow a toroid: start from each chain that is a toroid by itself, add one more
                        // chain orbit per round, keep the `beam` unions of highest genus, repeat.
                        // A union is recorded when it has more holes than what it grew from, or a
                        // different convex hull.
                        struct St { mem: Vec<usize>, reps: Vec<(usize, Xf)>, set: HashSet<(usize, [i64; 3])>, genus: i64, hull: String }
                        let mut beam_v: Vec<St> = Vec::new();
                        for (a, &i) in keep.iter().enumerate() {
                            if b[i].cls != "raqt" { continue; }
                            if let Some((genus, hull)) = try_union(&b[i].reps, "", &mut reasons, &|_, _| false) {
                                beam_v.push(St { mem: vec![a], reps: b[i].reps.clone(), set: sets[a].clone(), genus, hull });
                            }
                        }
                        let mut visited: HashSet<Vec<usize>> = HashSet::new();
                        for round in 0..rounds {
                            let mut nextv: Vec<St> = Vec::new();
                            for st in &beam_v {
                                for c in 0..keep.len() {
                                    if st.mem.contains(&c) || sets[c].is_subset(&st.set) { continue; }
                                    let mut mem = st.mem.clone();
                                    mem.push(c);
                                    mem.sort();
                                    if !visited.insert(mem.clone()) { continue; }
                                    let mut reps = st.reps.clone();
                                    for (s, x) in &b[keep[c]].reps[1..] { if !st.set.contains(&pkey(*s, x.t)) { reps.push((*s, *x)); } }
                                    if reps.len() >= 64 { continue; }
                                    let chain = mem.iter().map(|&m| b[keep[m]].desc.as_str()).collect::<Vec<_>>().join(" || ");
                                    let (g0, h0) = (st.genus, st.hull.clone());
                                    if let Some((genus, hull)) = try_union(&reps, &chain, &mut reasons, &|g, h| g > g0 || h != h0) {
                                        let mut set = st.set.clone();
                                        set.extend(sets[c].iter().cloned());
                                        nextv.push(St { mem, reps, set, genus, hull });
                                    }
                                }
                            }
                            if nextv.is_empty() { break; }
                            nextv.sort_by(|x, y| y.genus.cmp(&x.genus).then(x.reps.len().cmp(&y.reps.len())));
                            nextv.truncate(beam);
                            eprintln!("placement {} round {} unions kept {} best genus {} ({} chain orbits) hull {}", bi, round + 1, nextv.len(), nextv[0].genus, nextv[0].mem.len(), nextv[0].hull.replace('\t', " "));
                            beam_v = nextv;
                        }
                        continue;
                    }
                    let np = keep.len() * keep.len().saturating_sub(1) / 2;
                    if np > maxpairs { *reasons.entry("placement-skipped").or_insert(0) += 1; *reasons.entry("pairs-skipped").or_insert(0) += np as u64; continue; }
                    for a in 0..keep.len() { for c in a + 1..keep.len() {
                        let (pa, pc) = (&b[keep[a]], &b[keep[c]]);
                        if sets[c].is_subset(&sets[a]) || sets[a].is_subset(&sets[c]) { *reasons.entry("nested").or_insert(0) += 1; continue; }
                        let mut reps = pa.reps.clone();
                        for (s, x) in &pc.reps[1..] { if !sets[a].contains(&pkey(*s, x.t)) { reps.push((*s, *x)); } }
                        let chain = format!("{} [{}] || {} [{}]", pa.desc, pa.cls, pc.desc, pc.cls);
                        try_union(&reps, &chain, &mut reasons, &|_, h| shared.lock().unwrap().hulls.as_mut().map_or(true, |k| k.insert(h.to_string())));
                    } }
                }
                let mut sh = shared.lock().unwrap();
                for (k, v) in reasons { *sh.reasons.entry(k).or_insert(0) += v; }
            });
        }
    });
    let sh = shared.lock().unwrap();
    let mut rs: Vec<_> = sh.reasons.iter().collect();
    rs.sort();
    println!("done {:.1}s unique {} reasons {:?}", t0.elapsed().as_secs_f64(), sh.n, rs);
}
