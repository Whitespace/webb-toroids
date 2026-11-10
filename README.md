# webb-toroids

Search for Webb-style Stewart toroids: polyhedra that are (R) regular-faced, (A) free of coplanar
adjacent faces, (Q) quasi-convex (every convex hull edge is an edge of the toroid) and (T) of genus
at least 1 without self-intersection, where the convex hull is equilateral but may have irregular
faces. Only octahedral and icosahedral symmetry is searched.

## Build

Needs a stable Rust toolchain (edition 2021, no dependencies).

    cargo build --release

The binary is `target/release/webb-toroids`. It reads `solids.txt` from the current directory
(override with `--solids PATH`).

The Python scripts need Python 3 with `numpy` and `scipy`.

## Run

    target/release/webb-toroids --depth 1 --any 0 --threads 16 --out runs/t1

| flag | meaning | default |
| --- | --- | --- |
| `--depth D` | most links in a chain between the two axis pieces | 2 |
| `--any K` | how many of those links may be arbitrary solids; the rest must be entered and left through antiparallel faces | 1 |
| `--min-group N` | 24 keeps octahedral and icosahedral, 60 icosahedral only | 24 |
| `--near N` | also keep results with up to N hull edges missing (turns off the pre-screen) | 0 |
| `--only a,b,c` | restrict the library to these solids (names as in `solids.txt`) | all |
| `--seeds a,b,c` | restrict which solids may be the seed | all |
| `--seed-range A:B` | run seeds A to B-1 only, to split a run | all |
| `--threads T` | worker threads | 4 |
| `--out DIR` | output directory | `out` |

Output: `results.tsv` (one row per distinct toroid), `results.jsonl` (representative pieces and
group generators) and `off/NNNNNN.off` with a `.pieces` file naming the solid behind each face.
Progress and reject counts go to stderr.

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

## Included runs

- `runs/t1`: one through-link, all solids. Complete.
- `runs/t2c`: two through-links on a hand-picked set of 31 solids. Complete. Run with `--depth 2 --any 0 --only`
  `pentagonal-rotunda,pentagonal-prism,pentagonal-antiprism,tridiminished-icosahedron,augmented-tridiminished-icosahedron,metabidiminished-icosahedron,pentagonal-pyramid,gyroelongated-pentagonal-pyramid,elongated-pentagonal-pyramid,pentagonal-cupola,icosidodecahedron,dodecahedron,icosahedron,decagonal-prism,decagonal-antiprism,triangular-prism,cube,octahedron,tetrahedron,square-pyramid,square-cupola,triangular-cupola,square-antiprism,cuboctahedron,hexagonal-prism,octagonal-prism,truncated-tetrahedron,gyrobifastigium,elongated-square-pyramid,bilunabirotunda,triangular-hebesphenorotunda`
- `runs/t2`: two through-links, all solids. Stopped after 90 of 1847 seeds.
- `runs/webb4`: Webb's four piece types only, used as the sanity test.
