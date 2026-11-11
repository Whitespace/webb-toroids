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

Directly:

    target/release/webb-toroids --depth 1 --any 0 --out runs/t1

| flag | meaning | default |
| --- | --- | --- |
| `--depth D` | most links in a chain between the two axis pieces | 2 |
| `--any K` | how many of those links may be arbitrary solids; the rest must be entered and left through antiparallel faces | 1 |
| `--min-group N` | 24 keeps octahedral and icosahedral, 60 icosahedral only | 24 |
| `--near N` | also keep results with up to N hull edges missing (turns off the pre-screen) | 0 |
| `--max-radius R` | drop chains with a piece centre further than R edge lengths from the seed axis or the centre (Webb's toroid: 6.2) | none |
| `--max-pieces N` | reject candidates with more than N pieces (Webb's toroid: 152) | none |
| `--only a,b,c` | restrict the library to these solids (names as in `solids.txt`) | all |
| `--seeds a,b,c` | restrict which solids may be the seed | all |
| `--seed-range A:B` | run seeds A to B-1 only, to split a run | all |
| `--threads T` | worker threads | all CPUs |
| `--resume` | continue from `checkpoint.txt` in the output directory if there is one; refuses if the settings differ | off |
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

## Included runs

- `runs/t1`: one through-link, all solids. Complete.
- `runs/t2c`: two through-links on a hand-picked set of 31 solids. Complete. Run with `--depth 2 --any 0 --only`
  `pentagonal-rotunda,pentagonal-prism,pentagonal-antiprism,tridiminished-icosahedron,augmented-tridiminished-icosahedron,metabidiminished-icosahedron,pentagonal-pyramid,gyroelongated-pentagonal-pyramid,elongated-pentagonal-pyramid,pentagonal-cupola,icosidodecahedron,dodecahedron,icosahedron,decagonal-prism,decagonal-antiprism,triangular-prism,cube,octahedron,tetrahedron,square-pyramid,square-cupola,triangular-cupola,square-antiprism,cuboctahedron,hexagonal-prism,octagonal-prism,truncated-tetrahedron,gyrobifastigium,elongated-square-pyramid,bilunabirotunda,triangular-hebesphenorotunda`
- `runs/t2`: two through-links, all solids. Stopped after 90 of 1847 seeds.
- `runs/d3a0-r6.2-p152`: three through-links, all solids, `--max-radius 6.2 --max-pieces 152`. Stopped at 1731 of
  1847 seeds; continue with `mise run search 3 0 --max-radius 6.2 --max-pieces 152`. Same 33 results as `runs/d2` so far.
- `runs/d2a2-c31`: two arbitrary (bent) links on the 31-solid set of `t2c`, `--max-radius 6.2 --max-pieces 152`.
  Complete, 79 results. `scripts/check.py` agrees on 78; it reports 4 missing hull edges on 000019, where it splits 2 of
  the 30 octagonal hull faces, which looks like a hull-merging artefact and is not resolved.
- `runs/webb4`: Webb's four piece types only, used as the sanity test.
