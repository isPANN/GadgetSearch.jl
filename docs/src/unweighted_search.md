# Four-pin unweighted lattice search

`search_unweighted_gadgets` searches a finite lattice window for a four-pin
gadget whose reduced alpha tensor differs from the target by one constant. The
four pins and their outward directions must form the crossing geometry checked
by `check_gadget_geometry`.

One SAT formulation chooses the occupied sites, ordered pins, and outward rays
together. Supply one exact `window_shape`, `atom_count`, and tensor `offset`.
All vertices have unit weight, and edges follow the lattice blockade geometry.
The triangular lattice is the default; pass `Square()` to select KSG.

The target graph supplies the reduced alpha tensor and the final verifier.
The replacement must satisfy `α̃(replacement) = α̃(target) + offset`, including
the same `-Inf` entries. Every returned candidate is
independently checked by `is_gadget_replacement`, connectivity, and
`check_gadget_geometry`.

The first pin is fixed at the canonical origin. `canonical_shift` positions the
window relative to that pin, and `first_ray` selects its outward lattice
direction. On the triangular lattice, direction 1 is `(1, 0)` in axial
coordinates. The window must include the origin.

```julia
using GadgetSearch, Graphs

target = SimpleGraph(6)  # CROSS+EDGE; vertices 5 and 6 are internal
for edge in ((1, 2), (1, 5), (2, 6), (3, 5), (4, 6))
    add_edge!(target, edge...)
end

gadget = search_unweighted_gadgets(
    target, [1, 2, 3, 4], Triangular();
    window_shape=(3, 4), atom_count=9, offset=1,
    seconds=30, seed=1, canonical_shift=(-1, 1), first_ray=1,
)
if gadget !== nothing
    @assert is_gadget_replacement(
        target, gadget.replacement_graph, [1, 2, 3, 4],
        gadget.boundary_vertices,
    ) == (true, 1.0)
end
```

For CROSS alone, use the two edges `(1, 3)` and `(2, 4)`. This search requests
a 23-atom replacement with offset +7 and allows more solving time:

```julia
cross = SimpleGraph(4)
add_edge!(cross, 1, 3)
add_edge!(cross, 2, 4)

cross_gadget = search_unweighted_gadgets(
    cross, [1, 2, 3, 4], Triangular();
    window_shape=(8, 6), atom_count=23, offset=7,
    seconds=120, seed=1, canonical_shift=(-6, 3), first_ray=2,
)
if cross_gadget !== nothing
    @assert is_gadget_replacement(
        cross, cross_gadget.replacement_graph, [1, 2, 3, 4],
        cross_gadget.boundary_vertices,
    ) == (true, 7.0)
end
```

The function returns a verifier-accepted `UnweightedGadget` or `nothing`.
`nothing` means the specified instance was unsatisfiable or the SAT-solving
budget expired; those outcomes are not distinguished. `seconds` defaults to
600 and limits the solving loop, excluding CNF construction and final
verification. Schedule different counts, offsets, and symmetry settings
externally to search a larger space. One successful instance does not certify
a globally minimum atom count.

The former `search_unweighted_gadget_joint` API is now named
`search_unweighted_gadgets`. The old enumeration keywords, checkpoint API, and
`UnweightedSearchResult` have been removed. Callers now use the exact-instance
keywords above and access the returned gadget directly.

For a runnable example that searches both CROSS+EDGE and CROSS, verifies each
result, and exports separate JSON embeddings and PNG plots, see
[`examples/unweighted_Rydberg_example.jl`](https://github.com/isPANN/GadgetSearch.jl/blob/main/examples/unweighted_Rydberg_example.jl).

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
point under these rules and budgets is not a global minimum certificate. Each
fixed-frame solve is capped by `max_sat_conflicts`; capped calls are counted in
`unresolved_sat_evaluations`. If the neighborhood is exhausted while any such
call remains unresolved, the termination reason is `:sat_unknown`, not
`:rewrite_fixed_point`.

```julia
optimized = optimize_unweighted_gadget(
    gadget, [1, 2, 3, 4];
    min_vertices=4,
    max_sat_evaluations=256,
    max_sat_conflicts=100_000,
    host_radius=1,
)
```
