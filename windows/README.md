# Window library

A *window* here is a face of the convex hull that is not a face of any solid: a convex, planar
polygon with unit edges, like the irregular 12-gons of Webb's toroid. The `windows` binary lists
every such polygon that a ring of solids can frame, and for each rim edge which solids can sit
there. The aim is to start a toroid search from the hull: pick windows first, then fill in.

## What is in here

| directory | contents | command |
| --- | --- | --- |
| `open/` | fully open windows with 7 or more sides: nothing but rim edges lies in the plane | `windows --min-n 7 --no-faces-in-plane --out windows/open` |
| `covered13-direct/` | 13 or more sides, solid faces allowed in the plane, neighbouring rim pieces glued directly | `windows --min-n 13 --no-via --out windows/covered13-direct` |
| `covered13/` | the same with hidden corner pieces allowed (47 cpu-minutes) | `windows --min-n 13 --out windows/covered13` |

Each directory holds

- `windows.jsonl`: one polygon per line. `n`, `turns` (exterior angle after each edge, degrees),
  `poly` (corners, counter-clockwise seen from outside, first edge along +x), `rot` and `mirror`
  (symmetry of the turn sequence), `outline` (same number for a polygon and its mirror image),
  `rims` (number of closed rim walks, each counted once per starting edge), `options[i]` (ids of
  the states a piece on edge i can take) and `example` (one rim checked for overlaps: the pieces
  with their placement as 9 rotation entries by rows then 3 translation entries, window plane
  z = 0, outside is +z; `hidden` marks corner pieces that do not carry a rim edge).
- `states.tsv`: a state is a solid touching the plane along one directed edge at one tilt. `xf`
  places the solid with that edge from (0,0,0) to (1,0,0) and the plane at z = 0.
- `corners.tsv`: which state can follow which, the turn, and the hidden piece if there is one
  (`via_xf` places it in the frame of the `from` state).

## Query

    python scripts/windows.py windows/open list                # table of all polygons
    python scripts/windows.py windows/open list --n 12 --solid pentagonal-rotunda
    python scripts/windows.py windows/open show 27              # per-edge options of one window
    python scripts/windows.py windows/open off 27 /tmp/w.off    # example rim as OFF
    python scripts/check_windows.py windows/open                # independent check (numpy, scipy)

`mise run windows -- --min-n 12 --only a,b,c --out DIR` runs a search; flags: `--min-n`,
`--max-n`, `--only`, `--max-faces`, `--solids FILE`, `--no-faces-in-plane`, `--no-via`,
`--allow-coplanar`, `--threads`, `--example-budget`, `--full-tables` (write every usable state
and corner, not only those that occur in a window).

## Model

A rim piece lies under the plane and touches it along one edge (or along a whole face, when
faces in the plane are allowed; the rim then follows that face). Two neighbouring rim pieces
meet at a corner of the polygon and are either

1. glued face to face, on a face that contains the corner, or
2. both glued to a third, hidden piece, on two of its faces that contain the corner; the hidden
   piece touches the plane at the corner only.

The plane through the two rim edges must support all pieces of the corner. That fixes the tilt
of both rim pieces, so the corners form a finite graph on (solid, edge, tilt) states, and a
window is a closed walk whose turns add to 360 degrees and whose unit steps close up. The search
enumerates turn sequences (smallest turn first, so n at least N forces a first turn of at most
360/N) and keeps the set of walks that fit; it is exhaustive for the model. Afterwards each
polygon gets its per-edge options, pieces whose in-plane face would stick out of the polygon are
removed, and one overlap-free rim is searched as an example.

Not covered: rim pieces that only touch at the corner and are held by pieces further away, two
or more hidden pieces at one corner, corners of exactly 180 degrees, and (A) for hidden corners.
A window being in the library says nothing yet about a whole toroid existing around it.

## Results (2026-11-19, all 127 solids)

- Fully open windows have 8, 10 or 12 sides, never more: 17 + 9 + 2 polygons (12 + 6 + 2 up to
  mirror image). Allowing hidden corner pieces, or adding 7-, 9-, 11- and 12-gonal prisms and
  antiprisms, does not produce a 13-gon or larger.
- The two 12-gons are Webb's (turns 4 x 26.565 and 8 x 31.717, rim: rotunda, antiprism, prism,
  J63, ...) and one of cubic type (4 x 19.471 and 8 x 35.264).
- With pentagonal solids the full list is: Webb's octagon (the prism-free shell), Webb's 12-gon,
  a 10-gon half way between them (struts with a prism on one side only, chiral pair), a regular
  decagon, and one more octagon whose only rims overlap.
- From 13 sides on every window has solid faces lying in its plane (3 to 9 of them), and all are
  of cubic type (turns 10.53, 19.47, 35.26, 54.74, 70.53 degrees: square cupolas, octahedra,
  augmented triangular prisms): 66 polygons with 14, 15, 16, 18 and 24 sides when neighbouring
  rim pieces are glued directly.
- With hidden corner pieces allowed (`covered13/`): 325 polygons with 14, 15, 16, 18, 20, 21, 22
  and 24 sides, again all with solid faces in the plane; 310 have a checked example rim. The 15
  without one include regular 16- and 20-gons and other 18 and 22.5 degree polygons, the only
  ones that list pentagonal solids; treat those as unconfirmed.
- Stretching Webb's 12-gon to 14 sides (`scripts/bridge_chains.py`): the antiprism end of the
  rotunda arc accepts only a pentagonal prism (turn 31.72) or an icosahedron fragment (63.43).
  A 14-gon needs four pieces between the two antiprisms turning 126.87 degrees in total; of the
  20955 four-piece chains that fit, the flattest turns 180.
