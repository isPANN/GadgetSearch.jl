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
