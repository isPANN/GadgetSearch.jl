# Four-pin unweighted lattice search

`search_unweighted_gadgets` searches a finite lattice window for a four-pin
gadget whose reduced alpha tensor differs from the target by one constant. The
four pins and their outward directions must form the crossing geometry checked
by `check_crossing_frame`.

The target graph is used only to obtain the target reduced alpha tensor and for
the unchanged final verifier. The search chooses the concrete pin locations,
pin directions, and lattice sites. `window_side` controls the finite search
window. The triangular lattice is the default; pass `Square()` to select KSG.

The direct search enumerates the finite set of ordered pin locations and
outward lattice directions that satisfy the crossing geometry. For each layout,
atom count, and constant offset, a SAT instance chooses all occupied sites at
once. Its constraints enforce the target reduced alpha tensor and connectivity.
This does not require a known logical graph, a 28-site seed, or a historical
rewrite path.

`min_vertices:max_vertices` is examined in increasing order, so reaching a
larger atom count means all scheduled smaller cases have already been rejected.
`max_evaluations` bounds SAT solver calls, not graph mutations. For each fixed
layout, atom count, and offset, satisfying site selections are blocked and the
solver is called again until that case is unsatisfiable. Every returned
candidate is independently checked by both `is_gadget_replacement` and
`check_crossing_frame`.

`max_frame_evaluations` separately bounds pin-and-ray geometry checks. Increase
it together with `window_side` for exhaustive searches in larger windows.

### Checkpoint and resume

Pass `checkpoint_path` to persist deterministic frame, label-order, and offset
cursors after every `checkpoint_interval` SAT evaluations. Repeating the same
search call with the same path resumes from those cursors and restores its
evaluation and staged-filter counters. The target, lattice, vertex range, and
window size must match the saved search.

`read_unweighted_search_checkpoint(path)` returns the saved `lattice`,
`atom_count`, `window_shape`, `frame_cursor`, `order_cursor`, `offset_cursor`,
`evaluated`, `frame_evaluated`, `frame_candidates`, `first_lower_rejected`,
`second_lower_rejected`, and `full_solves` fields for progress reporting.

```julia
using GadgetSearch, Graphs

target = SimpleGraph(4)
add_edge!(target, 1, 3)
add_edge!(target, 2, 4)

result = search_unweighted_gadgets(
    target, [1, 2, 3, 4], Triangular();
    min_vertices=4,
    max_vertices=23,
    max_evaluations=100_000,
    max_frame_evaluations=10_000_000,
    max_results=1,
    window_side=8,
)
```

The keyword defaults are `min_vertices=5`, `max_vertices=13`,
`max_evaluations=2_000`, `max_frame_evaluations=1_000_000`, `max_results=1`,
and `window_side=4` for a four-pin target. The result contains the target,
ordered boundary, lattice name, accepted gadgets, number of SAT solver calls, and
one of four termination reasons: `:solution`, `:budget`, `:frame_budget`, or
`:search_space_exhausted`.

## Joint occupancy-and-frame search

`search_unweighted_gadget_joint` is the direct formulation used for blank-window
discovery. One CNF chooses the occupied sites, pins, and rays together. The
caller supplies one exact window shape, atom count, and tensor offset; the
function returns either a verifier-accepted `UnweightedGadget` or `nothing`.
`canonical_shift` and `first_ray` identify the symmetry-fixed first port for
this exact CNF instance. A large search should schedule those independent
instances externally with separate solver budgets.

```julia
gadget = search_unweighted_gadget_joint(
    target, [1, 2, 3, 4], Triangular();
    window_shape=(8, 8), atom_count=28, offset=10,
    seconds=600, seed=1, canonical_shift=(0, 0), first_ray=1,
)
```

## Rewrite optimization

`optimize_unweighted_gadget` is a local downstream atom-count optimizer. It uses
semantic rewrite rules rather than arbitrary vertex deletion:

- contract a two-edge boundary tail while promoting its endpoint to the pin;
- contract opposite leaf pins (`P1`–`P3` or `P2`–`P4`) as one paired rewrite;
- move one frame pin by one lattice step and let fixed-frame SAT re-synthesize
  every interior atom at a smaller atom count.

Every accepted step passes the unchanged reduced-alpha verifier and all four
crossing-frame checks. Direct rules are explored to a closure rather than
greedily committing to the first smaller graph. Fixed-frame SAT is then tried
from the direct descendants, so a smaller direct dead end does not hide a
rewrite-and-resynthesize path. The result includes the best gadget reached in
the selected one-step frame neighborhood, a replayable before/after rewrite
trace, the number of fixed-frame SAT calls, and its termination reason. A fixed
point under these rules and budgets is not a global minimum certificate.

```julia
optimized = optimize_unweighted_gadget(
    gadget, [1, 2, 3, 4];
    min_vertices=17,
    max_sat_evaluations=256,
    host_radius=1,
)
```
