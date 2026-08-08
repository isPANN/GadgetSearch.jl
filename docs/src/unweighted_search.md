# Unweighted gadget search

The unweighted search dynamically constructs concrete induced lattice patches.
Call it with `Square()` for KSG or `Triangular()` for the triangular lattice. It
does not enumerate subsets of a fixed rectangular canvas, and it never reports
an abstract graph without an embedding.

## Search state

A state consists of:

- a finite set of integer lattice coordinates;
- an ordered list of pin coordinates;
- one outward lattice ray for each pin.

The induced graph and boundary vertex indices are derived from those coordinates.
The boundary order is part of the state, so two patches with the same pins in a
different logical order are evaluated separately.

## Search loop

For four-pin targets, `search_unweighted_gadgets` starts from legal four-arm
frames whose arm lengths vary independently. It repeats four steps until it
exhausts the evaluation budget:

1. Propose geometric edits: add, remove, or relocate a site; move or swap pins;
   change a pin ray; extend an arm by two sites; or locally split a crowded
   non-pin site.
2. Add fresh legal frames so the search does not depend on one lineage.
3. Rebuild the induced lattice graph and compute its reduced alpha tensor.
4. Keep the best-scoring states plus a random exploration fraction for the next
   beam.

The coordinate plane is unbounded. The patch grows only where an action adds a
site, so increasing the allowed vertex count does not create a rectangular
combinatorial search space.

Every four-pin mutation must keep G1-G4 true before its tensor is evaluated.
Pins and rays therefore do not contribute a large post-hoc combination search.

The ranking score is lexicographic:

1. number of positions where only one tensor is infinite;
2. spread of the finite entry-wise offsets;
3. number of failed crossing-frame conditions G1-G4;
4. vertex and edge counts.

The first two terms guide candidates toward the verifier contract. For a four-pin
search, a result is returned only when all four geometric conditions also pass:
the interfaces are strict convex-hull vertices (G1), the two channels alternate
around the hull (G2), every ray points outward (G3), and the infinite exterior
corridors are empty, touch only their own pin, and are pairwise non-adjacent
(G4). The last two terms prefer smaller candidates when the earlier terms tie.
The reduced-alpha score never replaces `is_diff_by_constant`.

The default budget is 2,000 distinct tensor evaluations. Increase it only in a
controlled compute environment. Accepted candidates do not stop the run early;
the search keeps the best `max_results` complete crossing frames.

## Search trace

The returned `UnweightedSearchResult.trace` contains one
`UnweightedSearchRecord` per distinct tensor evaluation. Each record stores:

- the lattice type, occupied coordinates, ordered pin coordinates and rays,
  graph6 state, and derived boundary indices;
- the target graph6 state and target boundary in JSONL exports;
- its parent state and graph edit action;
- tensor-distance components;
- whether the state survived beam selection;
- whether the final verifier accepted it and, if so, the constant offset.

`save_unweighted_trace(path, result)` writes the same records as JSON Lines. This
keeps large future runs streamable and makes the trajectory directly consumable
from Python without serializing Julia graph objects.
