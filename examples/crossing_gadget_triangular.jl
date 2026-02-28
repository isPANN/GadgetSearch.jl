# # Search for Crossing Gadgets on Triangular Lattice UDG
#
# This example searches for **unweighted MIS gadget replacements** of the
# CROSS graph on triangular-lattice unit disk graphs (UDGs).
#
# Background (Theorem 3.7 in the paper):
#   Two graphs R and R' are valid MIS replacements of each other iff their
#   **reduced alpha tensors** α̃(R) and α̃(R') differ only by a constant offset.
#
# Pipeline:
#   Step 1 – Define the CROSS pattern graph (4 vertices, edges 1-3 and 2-4)
#   Step 2 – Load pre-generated triangular UDG datasets (pinset = [1,2,3,4])
#   Step 3 – Run search_unweighted_gadgets over each dataset
#   Step 4 – Report and visualise every found replacement

using GadgetSearch
using Graphs

# ============================================================
# Step 1: Define the CROSS pattern graph
# ============================================================
#
#   1       2
#    \     /
#     (cross)
#    /     \
#   3       4
#
# Edges: 1-3 and 2-4  (two paths that "cross" without sharing a middle vertex)

cross = SimpleGraph(4)
add_edge!(cross, 1, 3)
add_edge!(cross, 2, 4)
cross_boundary = [1, 2, 3, 4]

@info "CROSS graph: $(nv(cross)) vertices, $(ne(cross)) edges, boundary = $cross_boundary"

# ============================================================
# Step 2: Load triangular UDG datasets
# ============================================================
#
# Pre-generated files live in data/grid_udgs/.
# Naming convention: m{nx}n{ny}pad1_min{v_min}max{v_max}_direct4.g6
#   - nx×ny  inner grid size
#   - 4 boundary pins always at vertex indices [1,2,3,4]
#
# We search over increasing grid sizes and stop as soon as we have results.

datasets = [
    pkgdir(GadgetSearch, "data", "grid_udgs", "m3n3pad1_min3max9_direct4.g6"),
    pkgdir(GadgetSearch, "data", "grid_udgs", "m3n4pad1_min3max12_direct4.g6"),
    pkgdir(GadgetSearch, "data", "grid_udgs", "m3n5pad1_min3max15_direct4.g6"),
    pkgdir(GadgetSearch, "data", "grid_udgs", "m4n4pad1_min4max16_direct4.g6"),
]

all_results = UnweightedGadget[]

for path in datasets
    isfile(path) || (@warn "Dataset not found, skipping: $path"; continue)

    loader = GraphLoader(path; pinset=cross_boundary)
    @info "Searching $(basename(path))  ($(length(loader)) graphs)"

    t = @elapsed results = search_unweighted_gadgets(cross, cross_boundary, loader)

    @info "  → found $(length(results)) replacement(s) in $(round(t; digits=2))s"
    append!(all_results, results)
end

# ============================================================
# Step 3: Report results
# ============================================================

println()
println("=" ^ 60)
println("  CROSSING GADGET SEARCH  –  SUMMARY")
println("=" ^ 60)
println("Pattern:  CROSS graph  (4 vertices, edges 1-3 and 2-4)")
println("Space:    triangular-lattice UDGs, pinset = [1,2,3,4]")
println("Results:  $(length(all_results)) replacement gadget(s) found")
println()

for (i, r) in enumerate(all_results)
    println("── Gadget #$i ──────────────────────────────────────────")
    println("  vertices       : $(nv(r.replacement_graph))")
    println("  edges          : $(ne(r.replacement_graph))")
    println("  boundary       : $(r.boundary_vertices)")
    println("  constant offset: $(r.constant_offset)")
    if !isnothing(r.pos)
        println("  positions      : $(r.pos)")
    end
    println()
end

# ============================================================
# Step 4: Verify correctness and visualise
# ============================================================

out_dir = pkgdir(GadgetSearch, "examples")

for (i, r) in enumerate(all_results)
    # --- correctness check ---
    valid, offset = is_gadget_replacement(
        cross, r.replacement_graph,
        cross_boundary, r.boundary_vertices
    )
    status = valid ? "✓ valid" : "✗ INVALID"
    @info "Gadget #$i verification: $status  (offset = $offset)"

    # --- visualise replacement graph ---
    pos = isnothing(r.pos) ? nothing : r.pos
    out_path = joinpath(out_dir, "crossing_gadget_$(i).svg")
    GadgetSearch.plot_graph(r.replacement_graph, out_path;
        pos=pos,
        plot_size=400,
        margin=40,
        vertex_size=12,
        vertex_label_size=14,
        edge_width=2
    )
end

