using GadgetSearch
using Graphs
using GraphIO
using Test

# Build CROSS graph: 4 vertices, edges 1-3 and 2-4
function _cross_graph()
    g = SimpleGraph(4)
    add_edge!(g, 1, 3); add_edge!(g, 2, 4)
    return g
end

# Build BATOIDEA graph: 11-vertex replacement for CROSS (Figure 6 in the paper)
function _batoidea_graph()
    g = SimpleGraph(11)
    add_edge!(g, 1, 5);  add_edge!(g, 1, 9)
    add_edge!(g, 2, 5);  add_edge!(g, 2, 6);  add_edge!(g, 2, 7)
    add_edge!(g, 3, 8)
    add_edge!(g, 4, 9);  add_edge!(g, 4, 10); add_edge!(g, 4, 11)
    add_edge!(g, 5, 6);  add_edge!(g, 5, 9);  add_edge!(g, 5, 10)
    add_edge!(g, 6, 7);  add_edge!(g, 6, 9);  add_edge!(g, 6, 10); add_edge!(g, 6, 11)
    add_edge!(g, 7, 8);  add_edge!(g, 7, 10); add_edge!(g, 7, 11)
    add_edge!(g, 8, 11)
    add_edge!(g, 9, 10)
    add_edge!(g, 10, 11)
    return g
end

# Encode a graph to a Graph6 string (strips the >>graph6<< header)
function _to_g6(g)
    return GraphIO.Graph6._graphToG6String(g)[11:end]
end

@testset "make_unweighted_filter: CROSS → BATOIDEA" begin
    cross = _cross_graph()
    batoidea = _batoidea_graph()

    # Build a loader containing only BATOIDEA
    loader = GraphLoader(GraphDataset([_to_g6(batoidea)]))

    filter_fn = make_unweighted_filter(cross, [1, 2, 3, 4])
    result = filter_fn(batoidea, nothing, nothing)

    @test result !== nothing
    @test result isa UnweightedGadget
    @test result.boundary_vertices == [1, 2, 3, 4]
    @test result.constant_offset == 2.0
end

@testset "search_unweighted_gadgets: finds BATOIDEA as replacement for CROSS" begin
    cross = _cross_graph()
    batoidea = _batoidea_graph()

    # Loader with two graphs: CROSS itself (offset 0) and BATOIDEA (offset +2)
    loader = GraphLoader(GraphDataset([_to_g6(cross), _to_g6(batoidea)]),
                         pinset=[1, 2, 3, 4])

    results = search_unweighted_gadgets(cross, [1, 2, 3, 4], loader)

    @test length(results) >= 1
    # BATOIDEA should appear with constant_offset == 2.0
    @test any(r -> r.constant_offset == 2.0, results)
end

# ---------------------------------------------------------------------------
# Self-replacement
# ---------------------------------------------------------------------------

@testset "make_unweighted_filter: self-replacement has offset 0" begin
    cross = _cross_graph()
    filter_fn = make_unweighted_filter(cross, [1, 2, 3, 4])

    result = filter_fn(cross, nothing, [1, 2, 3, 4])

    @test result !== nothing
    @test result isa UnweightedGadget
    @test result.constant_offset == 0.0
    @test result.boundary_vertices == [1, 2, 3, 4]
end

# ---------------------------------------------------------------------------
# Early-exit: candidate has too few vertices
# ---------------------------------------------------------------------------

@testset "make_unweighted_filter: candidate with too few vertices returns nothing" begin
    cross = _cross_graph()
    filter_fn = make_unweighted_filter(cross, [1, 2, 3, 4])  # k = 4

    # A 3-vertex triangle has nv < k=4, so the filter must exit immediately
    triangle = SimpleGraph(3)
    add_edge!(triangle, 1, 2); add_edge!(triangle, 2, 3); add_edge!(triangle, 1, 3)

    @test filter_fn(triangle, nothing, nothing) === nothing
end

# ---------------------------------------------------------------------------
# Keyword: limit
# ---------------------------------------------------------------------------

@testset "search_unweighted_gadgets: limit keyword caps examined graphs" begin
    cross    = _cross_graph()
    batoidea = _batoidea_graph()

    # All three graphs are valid replacements for CROSS (offsets 0, 2, 0)
    loader = GraphLoader(
        GraphDataset([_to_g6(cross), _to_g6(batoidea), _to_g6(cross)]),
        pinset=[1, 2, 3, 4]
    )

    results_all = search_unweighted_gadgets(cross, [1, 2, 3, 4], loader)
    results_lim = search_unweighted_gadgets(cross, [1, 2, 3, 4], loader; limit=1)

    @test length(results_all) == 3   # all three examined and matched
    @test length(results_lim) == 1   # only the first graph was examined
end

# ---------------------------------------------------------------------------
# Keyword: max_results
# ---------------------------------------------------------------------------

@testset "search_unweighted_gadgets: max_results stops after finding N results" begin
    cross = _cross_graph()

    # Three identical CROSS graphs — each is a valid replacement (offset 0)
    loader = GraphLoader(
        GraphDataset([_to_g6(cross), _to_g6(cross), _to_g6(cross)]),
        pinset=[1, 2, 3, 4]
    )

    results = search_unweighted_gadgets(cross, [1, 2, 3, 4], loader; max_results=2)

    @test length(results) == 2
end

# ---------------------------------------------------------------------------
# Integration: triangular lattice UDG dataset
# ---------------------------------------------------------------------------

@testset "search_unweighted_gadgets on triangular lattice UDG dataset" begin
    cross = _cross_graph()

    # Use the pre-generated 3×3 triangular UDG dataset (boundary pins at [1,2,3,4])
    data_path = pkgdir(GadgetSearch, "data", "grid_udgs", "m3n3pad1_min3max9_direct4.g6")
    if !isfile(data_path)
        @warn "Dataset not found, skipping integration test: $data_path"
        return
    end
    loader = GraphLoader(data_path; pinset=[1, 2, 3, 4])

    results = search_unweighted_gadgets(cross, [1, 2, 3, 4], loader)

    # --- structural checks ---
    @test all(r -> r isa UnweightedGadget,              results)
    @test all(r -> length(r.boundary_vertices) == 4,    results)
    @test all(r -> nv(r.replacement_graph) >= 4,        results)
    # replacement graph must have the declared boundary vertices as valid indices
    @test all(r -> maximum(r.boundary_vertices) <= nv(r.replacement_graph), results)

    # --- correctness: every result must pass is_gadget_replacement ---
    for r in results
        valid, offset = is_gadget_replacement(
            cross, r.replacement_graph,
            [1, 2, 3, 4], r.boundary_vertices
        )
        @test valid
        @test offset ≈ r.constant_offset
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Multi-target search tests
# ═══════════════════════════════════════════════════════════════════════════════

# Helper: three CROSS variants (all perfect matchings on 4 vertices)
function _cross_b()
    g = SimpleGraph(4); add_edge!(g, 1, 2); add_edge!(g, 3, 4); return g
end
function _cross_c()
    g = SimpleGraph(4); add_edge!(g, 1, 4); add_edge!(g, 2, 3); return g
end

# ---------------------------------------------------------------------------
# inf_mask
# ---------------------------------------------------------------------------

@testset "inf_mask" begin
    # All finite → mask = 0
    @test inf_mask([1.0, 2.0, 3.0, 4.0]) == UInt(0)

    # First element -Inf → bit 0 set
    @test inf_mask([-Inf, 2.0, 3.0, 4.0]) == UInt(1)

    # Third element -Inf → bit 2 set
    @test inf_mask([1.0, 2.0, -Inf, 4.0]) == UInt(4)

    # Multiple -Inf entries
    @test inf_mask([-Inf, 2.0, -Inf, 4.0]) == UInt(5)  # bits 0 and 2

    # CROSS α̃: check that two isomorphic targets produce the same mask
    cross = _cross_graph()
    cross_reduced = GadgetSearch.content.(GadgetSearch.calculate_reduced_alpha_tensor(cross, [1,2,3,4]))
    m = inf_mask(cross_reduced)
    @test m != UInt(0)  # some entries must be -Inf
    # Self-consistency: mask bits match actual -Inf positions
    for (i, v) in enumerate(cross_reduced)
        @test isinf(v) == ((m >> (i-1)) & 1 == 1)
    end
end

# ---------------------------------------------------------------------------
# pins_prefilter
# ---------------------------------------------------------------------------

@testset "pins_prefilter" begin
    # Graph where all 4 pins are connected → pass
    cross = _cross_graph()
    @test pins_prefilter(cross, [1,2,3,4]) == false  # pins 1,3 and 2,4 are in separate components!

    batoidea = _batoidea_graph()
    @test pins_prefilter(batoidea, [1,2,3,4]) == true  # all connected

    # Graph with an isolated pin → fail
    g = SimpleGraph(5)
    add_edge!(g, 1, 5); add_edge!(g, 2, 5); add_edge!(g, 3, 5)
    # Pin 4 has degree 0
    @test pins_prefilter(g, [1,2,3,4]) == false

    # All pins isolated → fail
    g2 = SimpleGraph(4)
    @test pins_prefilter(g2, [1,2,3,4]) == false
end

# ---------------------------------------------------------------------------
# make_multi_target_filter: basic match
# ---------------------------------------------------------------------------

@testset "make_multi_target_filter: CROSS variants → BATOIDEA" begin
    cross = _cross_graph()
    cross_b = _cross_b()
    cross_c = _cross_c()
    batoidea = _batoidea_graph()

    targets = Tuple{SimpleGraph{Int}, Vector{Int}}[
        (cross,   [1,2,3,4]),
        (cross_b, [1,2,3,4]),
        (cross_c, [1,2,3,4]),
    ]

    filter_fn = make_multi_target_filter(targets)

    # BATOIDEA should match CROSS (target index 1) with offset 2.0
    result = filter_fn(batoidea, nothing, [1,2,3,4])
    @test result !== nothing
    @test result isa MultiTargetResult
    @test result.target_index == 1
    @test result.gadget.constant_offset == 2.0
end

# ---------------------------------------------------------------------------
# make_multi_target_filter: self-replacement
# ---------------------------------------------------------------------------

@testset "make_multi_target_filter: self-replacement (prefilter=false)" begin
    cross = _cross_graph()
    targets = Tuple{SimpleGraph{Int}, Vector{Int}}[(cross, [1,2,3,4])]

    # CROSS has 2 components ({1,3},{2,4}), so prefilter rejects it.
    # Use prefilter=false to test pure tensor matching.
    filter_fn = make_multi_target_filter(targets; prefilter=false)

    result = filter_fn(cross, nothing, [1,2,3,4])
    @test result !== nothing
    @test result.target_index == 1
    @test result.gadget.constant_offset == 0.0

    # With prefilter=true, CROSS itself is correctly rejected (disconnected pins)
    filter_fn2 = make_multi_target_filter(targets; prefilter=true)
    @test filter_fn2(cross, nothing, [1,2,3,4]) === nothing
end

# ---------------------------------------------------------------------------
# make_multi_target_filter: prefilter rejects disconnected pins
# ---------------------------------------------------------------------------

@testset "make_multi_target_filter: prefilter rejects disconnected" begin
    cross = _cross_graph()
    targets = Tuple{SimpleGraph{Int}, Vector{Int}}[(cross, [1,2,3,4])]

    filter_fn = make_multi_target_filter(targets; prefilter=true)

    # A graph where pins are in 2 separate components → rejected
    g = SimpleGraph(6)
    add_edge!(g, 1, 5); add_edge!(g, 3, 5)  # component {1,3,5}
    add_edge!(g, 2, 6); add_edge!(g, 4, 6)  # component {2,4,6}
    @test filter_fn(g, nothing, [1,2,3,4]) === nothing

    # Same graph, but prefilter disabled → goes through to tensor check
    filter_fn_no = make_multi_target_filter(targets; prefilter=false)
    result = filter_fn_no(g, nothing, [1,2,3,4])
    # Will compute tensor but won't match CROSS → nothing
    @test result === nothing
end

# ---------------------------------------------------------------------------
# search_multi_target_gadgets: integration
# ---------------------------------------------------------------------------

@testset "search_multi_target_gadgets: finds BATOIDEA" begin
    cross   = _cross_graph()
    cross_b = _cross_b()
    cross_c = _cross_c()
    batoidea = _batoidea_graph()

    targets = Tuple{SimpleGraph{Int}, Vector{Int}}[
        (cross,   [1,2,3,4]),
        (cross_b, [1,2,3,4]),
        (cross_c, [1,2,3,4]),
    ]

    loader = GraphLoader(
        GraphDataset([_to_g6(batoidea), _to_g6(cross)]),
        pinset=[1,2,3,4]
    )

    results = search_multi_target_gadgets(targets, loader)
    @test length(results) >= 1
    @test all(r -> r isa MultiTargetResult, results)

    # BATOIDEA should match target 1 (CROSS) with offset 2.0
    bat_results = filter(r -> r.gadget.constant_offset == 2.0, results)
    @test length(bat_results) >= 1
    @test bat_results[1].target_index == 1
end

# ---------------------------------------------------------------------------
# search_multi_target_gadgets: max_results
# ---------------------------------------------------------------------------

@testset "search_multi_target_gadgets: max_results" begin
    cross = _cross_graph()
    batoidea = _batoidea_graph()
    targets = Tuple{SimpleGraph{Int}, Vector{Int}}[(cross, [1,2,3,4])]

    # Use BATOIDEA (connected graph) so prefilter doesn't reject
    loader = GraphLoader(
        GraphDataset([_to_g6(batoidea), _to_g6(batoidea), _to_g6(batoidea)]),
        pinset=[1,2,3,4]
    )
    results = search_multi_target_gadgets(targets, loader; max_results=2)
    @test length(results) == 2
end

