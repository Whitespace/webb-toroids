# webb-toroids

Search for Webb-style Stewart toroids: polyhedra that are (R) regular-faced, (A) free of coplanar
adjacent faces, (Q) quasi-convex (every convex hull edge is an edge of the toroid) and (T) of genus
at least 1 without self-intersection, where the convex hull is equilateral but may have irregular
faces. Only octahedral and icosahedral symmetry is searched.

## Build

With [mise](https://mise.jdx.dev), which installs Rust and Python and creates `.venv`:

    mise install
    mise run build

Other tasks: `mise run check runs/t1/off/*.off` and `mise run viewer` (both install numpy and
scipy from `requirements.txt` first).

Without mise it needs a stable Rust toolchain (edition 2021, no crate dependencies):

    cargo build --release

The binary is `target/release/webb-toroids`. It reads `solids.txt` from the current directory
(override with `--solids PATH`).

`.cargo/config.toml` builds for the CPU of the build machine (`-C target-cpu=native`), so the
binary is not portable to older CPUs; delete that file for a generic build. The hot loops in
`src/simd.rs` are written as fixed-width lane blocks over structure-of-arrays vertex tables and
are vectorised by the compiler: NEON on Apple Silicon, AVX2 or AVX-512 on x86-64. There are no
intrinsics and no crate dependencies.

The Python scripts need Python 3 with `numpy` and `scipy`.

## Run

The resumable way, as a mise task:

    mise run search 2 0                  # depth 2, no arbitrary links -> runs/d2a0
    mise run search 3 0 --min-group 60   # extra flags pass through   -> runs/d3a0-g60

The task builds first, names the output directory from its arguments, uses every CPU, and
resumes: stop it with Ctrl-C and run the same command again to continue. Progress goes to the
terminal and to `run.log` in the output directory.

Add `--push` to commit the run's directory and push it when the run completes (`--commit` to
commit only; `SEARCH_PUSH=1` in the environment does the same as `--push`). A run that was
interrupted is not committed.

Directly:

    target/release/webb-toroids --depth 1 --any 0 --out runs/t1

| flag | meaning | default |
| --- | --- | --- |
| `--depth D` | most links in a chain between the two axis pieces | 2 |
| `--any K` | how many of those links may be arbitrary solids; the rest must be entered and left through antiparallel faces | 1 |
| `--min-group N` | 24 keeps octahedral and icosahedral, 60 icosahedral only | 24 |
| `--near N` | also keep results with up to N hull edges missing (turns off the pre-screen) | 0 |
| `--links a,b,c` | solids allowed as links; `struts` makes the free links prisms and antiprisms (cube and octahedron included) passed base to base, so chains are straight struts; `--any K` then allows K other pieces in the chain as bends. End pieces are not restricted | all |
| `--max-faces N` | use only solids with at most N faces | all |
| `--max-radius R` | drop chains with a piece centre further than R edge lengths from the seed axis or the centre (Webb's toroid: 6.2) | none |
| `--max-pieces N` | reject candidates with more than N pieces (Webb's toroid: 152) | none |
| `--only a,b,c` | restrict the library to these solids (names as in `solids.txt`) | all |
| `--seeds a,b,c` | restrict which solids may be the seed | all |
| `--seed-range A:B` | run seeds A to B-1 only, to split a run | all |
| `--threads T` | worker threads | all CPUs |
| `--resume` | continue from `checkpoint.txt` in the output directory if there is one; refuses if the settings differ | off |
| `--no-join` | visit the last chain level link by link, as before the join (slower; same results, for checking) | off |
| `--no-share` | test every hit on its own instead of sharing the chain-prefix screens between hits with the same group (slower; same results, for checking) | off |
| `--new-hulls FILE` | record one toroid per convex hull, and none on a hull listed in FILE (lines of group, circumradius, inradius, hull face sizes, as `scripts/known_hulls.py` prints them from `runs/`); the rest are counted as `known-hull`. `mise run search` and the queue worker pass this by default, `mise run search ... --all-fillings` turns it off | off |
| `--partials` | also write every collision-free hit, toroid or not, to `partials.tsv` (turns the hull screens off, so slower; large file) | off |
| `--combine DIR` | no search: union pairs of chains from `DIR/partials.tsv` that share seed, group and centre, validate each, write `DIR/combined/` (`--max-pairs N` skips placements with more pairs, default 2,000,000) | |
| `--partial-keys FILE` | with `--partials`: write only placements listed in FILE (lines of seed index, axis index and the three key numbers, i.e. the first five columns of `partials.tsv`) | all |
| `--base-raqt`, `--rounds N`, `--beam K` | with `--combine`: grow toroids instead of pairing everything. Start from each chain that is a toroid by itself, add one chain orbit per round, keep the K unions of highest genus; a union is recorded when it has more holes than what it grew from or a different convex hull | off, 1, 8 |
| `--third L` | with `--combine` (pair mode): two chains from one seed give up to three kinds of node but only two kinds of strut. For every pair whose union is free of collisions, also run straight struts of up to L prisms and antiprisms out of every face of every piece and validate each | 0 (off) |
| `--frames` | also record toroids that pass the hull screens (every hull edge has unit length) but miss some hull edge; the screens stay on, unlike `--near` | off |
| `--decorate FILE --id N` | no search: grow result N of a `results.jsonl` outward. Level 1 puts every solid on every face of every representative piece and takes its orbit under the result's own group; each further level (up to `--depth`) adds another such decoration or a solid on a piece added before. Records only (R)(A)(Q)(T) toroids whose hull differs from the base's and from `--new-hulls`. A decoration need not end on an axis, so this reaches shapes the chain search cannot. Depth 2 on Webb's toroid takes 30 s on 4 cores; `--only` and `--max-faces` limit the solids added | |
| `--dec-beam K` | with `--decorate`: beam search instead of everything. Only decorations that reach beyond the old hull; the K states of a level with the smallest share of non-toroid hull edges are extended. Heuristic | 0 (off) |
| `--out DIR` | output directory | `out` |

Output: `results.tsv` (one row per distinct toroid), `results.jsonl` (representative pieces and
group generators), `off/NNNNNN.off` with a `.pieces` file naming the solid behind each face,
`run.log` (progress and reject counts, also on stderr) and `checkpoint.txt` (settings, completed
seeds and result signatures). Without `--resume` an existing output directory is overwritten.
A seed interrupted part-way is redone on resume; its results are not duplicated. After a resume
the reject counts on the `done` line cover the last session only. `WATCH=1` prints live counters
every 10 s.

Check a result independently (scipy convex hull):

    python3 scripts/check.py runs/t1/off/*.off

Rebuild the viewer page from chosen results:

    python3 viewer/build_viewer.py viewer/picks.json   # writes viewer/index.html

`solids.txt` holds the 5 Platonic, 13 Archimedean and 92 Johnson solids plus prisms and
antiprisms up to 10 sides, with mirror images of the chiral ones (127 in all), at unit edge.
`scripts/prep.py DATA_DIR solids.txt` regenerates it from per-solid JSON files (vertices, edges,
faces, as in tesseralis/polyhedra-viewer) and polishes the coordinates to about 1e-14.

## Background queue

Instead of running jobs by hand, a worker can chip away at a backlog on a couple of cores:

    brew install terminal-notifier   # optional, for notifications with a thumbnail
    mise run queue:install           # macOS: installs and starts a launchd agent, restarts at login
    mise run queue:status            # the queue, what is running, recent log
    mise run queue:add -- 4 1 --links struts --max-radius 6.2 --max-pieces 152
    mise run queue:uninstall         # stop it; the running job resumes next time
    mise run queue:reinstall         # uninstall then install: restart on the latest build
    mise run queue:push              # commit and push runs/ as it is now, the running job included
    mise run queue:run               # or run the worker in the foreground, any OS

`queue/jobs.txt` is the backlog, one `mise run search` argument list per line, worked top to
bottom. It is tracked in git and the worker pulls before each job, so jobs pushed to the repo
are picked up without touching the machine. A job is finished once its `run.log` has a `done`
line; finished runs are committed locally. Everything is resumable: stopping the worker, a
reboot or a crash loses at most the seeds in flight.

Settings go in `queue/config.local` (see `queue/config.example`): `THREADS` (default 2),
`BACKGROUND` (macOS background QoS through `taskpolicy -b`, which keeps the work on the
efficiency cores), `PULL`, `PUSH` (push finished runs, off by default), `NOTIFY`.

`mise run queue:notify-test` sends one notification with a thumbnail straight away and says which
backend it used, so you can check the setup without waiting for a find.

While a job runs, `scripts/queue_notify.py` looks once a minute for results that no other run
under `runs/` has, and posts one desktop notification for them with a thumbnail of the simplest
(`scripts/render.py`, a small numpy rasteriser). Thumbnails and logs are in `queue/state/`.

Not tested on macOS by the author of this section: the launchd agent, `taskpolicy` and the
notification itself. The worker loop, naming, novelty check, thumbnails and commits were
tested on Linux. If the repository lives under Documents, Desktop or Downloads, macOS may
block a launchd agent from reading it until bash is given Full Disk Access.

## Profiling

    mise run profile            # bent chains, 2 threads, 20 s of samples
    mise run profile struts 2 60

builds a second binary under `target/profiling` (line tables on, feature `prof` keeps the
mid-level functions out of line so they show up separately), runs a fixed workload, samples it
(`sample` on macOS, `perf` on Linux) and writes a text report to `profiles/`, with the search's
own per-reason cpu counters at the end. Workloads: `bent`, `struts`, `straight`.

## How it works

A seed solid sits on one of its own rotation axes of order 3, 4 or 5. A chain of solids glued
face to face leaves one of its faces. The chain is a hit when its last solid has one of its own
axes meeting the seed's axis at an angle that occurs between axes of the octahedral or
icosahedral group. The two axis rotations generate the group and the orbit of the chain is the
candidate. Candidates are then checked for overlap (separating axes, f64 with 1e-7 tolerance),
improper contact, closed manifoldness, genus, (A) and (Q).

Almost all the time goes into rejecting candidates, so the cheap necessary conditions for (Q)
run first and never build the candidate: the orbit is kept as piece centres, and each test
turns its directions into a solid's own frame and runs on that solid's vertex table. In order:
inscribed spheres overlapping; the cone of unit edges at the outermost vertex; then, edge by
edge from the outside in, the hull faces through that edge must have unit sides. Edges that a
piece's own symmetry carries onto one another are tested once. Only survivors are placed,
collision-tested and meshed.

The last chain level is a join rather than a loop. Links that leave through the same polygon
with the same rotation see the seed axis in the same direction, so they are grouped; the end
pieces' axis directions are stored once each. For each group one vector pass bins the angle
between the seed axis and every stored direction, and only pairs at an allowed angle whose axes
are also coplanar are looked at further. The pairs that survive are exactly the ones the plain
loop would have reached, in the same order (`--no-join` runs the plain loop for comparison).
Directions that coincide across groups, or are opposite, are merged into one table per polygon,
so the vector pass runs once per chain node over about 2.5 times fewer entries.

Tables of 8192 directions or more are stored cell by cell of a cube-map grid, about 32 to a cell.
The vector pass then runs over the cell centres against the allowed angle bins widened by the
cell size, and only cells that can hold a match are opened, on a compact single-precision copy.
Seeds that are all scan (many have no hit at all on all-solid runs) take a third to a quarter of
the time. `NO_GRID=1` gives the full pass for comparison; results and reject counts are the same.

A seed's nearest copy sits on the neighbouring axis of the same order, so below a closed-form
distance from the centre their inscribed spheres overlap whatever the chain does. Those centres
are rejected before the group is built (`seed-overlap`), and a seed for which that distance
exceeds `--max-radius` is skipped (`seed-dead`).

Many hits at one chain node generate the same group, and most of them fail on what the chain
so far does under that group, whatever is attached afterwards. So per node and group the
prefix is tested once: if its orbit overlaps itself or is too big, every hit with that group is
dropped; if it fails the outermost-vertex screen, a hit is dropped when its remaining pieces lie
closer to the centre than that vertex and put no vertex one edge length from it (then the full
orbit fails the same screen at the same vertex). Everything else takes the normal path, so the
results are the same (`--no-share` to compare; `SHARE_CHECK=1` asserts that hits sharing a cache
entry really have the same group). In the `done` line these show as `prefix-overlap` and
`prefix-q`.

The same cache keeps the prefix's orbit. A hit that is not dropped outright runs the screens on
that orbit plus its own remaining pieces, which is most of the per-hit cost saved; one that passes
them is validated again from scratch, so nothing recorded depends on the cache. One visible side
effect: when a link's centre lies on a rotation axis, the cached and the fresh orbit keep
differently turned copies of it, and a few rejections per 10,000 move between `q-prescreen` and
`q-proxy`. Both are rejections either way.

## Included runs

- `runs/t1`: one through-link, all solids. Complete.
- `runs/t2c`: two through-links on a hand-picked set of 31 solids. Complete. Run with `--depth 2 --any 0 --only`
  `pentagonal-rotunda,pentagonal-prism,pentagonal-antiprism,tridiminished-icosahedron,augmented-tridiminished-icosahedron,metabidiminished-icosahedron,pentagonal-pyramid,gyroelongated-pentagonal-pyramid,elongated-pentagonal-pyramid,pentagonal-cupola,icosidodecahedron,dodecahedron,icosahedron,decagonal-prism,decagonal-antiprism,triangular-prism,cube,octahedron,tetrahedron,square-pyramid,square-cupola,triangular-cupola,square-antiprism,cuboctahedron,hexagonal-prism,octagonal-prism,truncated-tetrahedron,gyrobifastigium,elongated-square-pyramid,bilunabirotunda,triangular-hebesphenorotunda`
- `runs/t2`: two through-links, all solids. Stopped after 90 of 1847 seeds.
- `runs/d2a2-c31`: two arbitrary (bent) links on the 31-solid set of `t2c`, `--max-radius 6.2 --max-pieces 152`.
  Complete, 79 results. `scripts/check.py` agrees on 78; it reports 4 missing hull edges on 000019, where it splits 2 of
  the 30 octagonal hull faces, which looks like a hull-merging artefact and is not resolved.
- `runs/d2a2-f14-r6.2-p152`: two bent links on the 47 solids with at most 14 faces (`mise run search 2 2 --max-faces 14
  --max-radius 6.2 --max-pieces 152`). Complete, 120 results, all pass `scripts/check.py`.
- `runs/d3a0-r6.2-p152`: complete (finished on a second machine), 33 results, the same as `runs/d2`.
- `runs/d8a0-struts-r6.2-p152`: straight struts of up to 8 prisms and antiprisms between any two end pieces
  (`mise run search 8 0 --links struts --max-radius 6.2 --max-pieces 152`). Complete in seconds, 14 results, all already in
  `runs/d2`. Struts up to 12 long with radius up to 12 and no piece limit give the same 14.
- `runs/d3a1-struts-r6.2-p152`: straight struts with one bend piece. Complete, 148 results, 66 not in any other
  run, none with a new hull. `scripts/check.py` agrees on 147 and disputes 000038 (the same shape it disputes in `d2a2-c31`).
- `runs/webb4`: Webb's four piece types only, used as the sanity test.

## Window library

`windows/` holds a library of hull windows (convex unit-edged polygons such as the 12-gons of
Webb's toroid) with the solids that can frame them, built by `src/bin/windows.rs`. See
`windows/README.md`.
