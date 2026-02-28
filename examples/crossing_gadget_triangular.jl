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
#   Step 2 – Verify a known replacement (BATOIDEA) using make_unweighted_filter
#   Step 3 – Search a small triangular UDG dataset with search_unweighted_gadgets
#   Step 4 – (offline) Full search over all pre-generated datasets

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

# ============================================================
# Step 2: Verify a known replacement using make_unweighted_filter
# ============================================================
#
# The BATOIDEA graph (11 vertices, Figure 6 in the paper) is a known
# crossing gadget replacement.  We pre-compute the filter once and
# then check BATOIDEA – this should return a valid UnweightedGadget
# with constant_offset == 2.0.

batoidea = SimpleGraph(11)
add_edge!(batoidea, 1, 5);  add_edge!(batoidea, 1, 9)
add_edge!(batoidea, 2, 5);  add_edge!(batoidea, 2, 6);  add_edge!(batoidea, 2, 7)
add_edge!(batoidea, 3, 8)
add_edge!(batoidea, 4, 9);  add_edge!(batoidea, 4, 10); add_edge!(batoidea, 4, 11)
add_edge!(batoidea, 5, 6);  add_edge!(batoidea, 5, 9);  add_edge!(batoidea, 5, 10)
add_edge!(batoidea, 6, 7);  add_edge!(batoidea, 6, 9);  add_edge!(batoidea, 6, 10); add_edge!(batoidea, 6, 11)
add_edge!(batoidea, 7, 8);  add_edge!(batoidea, 7, 10); add_edge!(batoidea, 7, 11)
add_edge!(batoidea, 8, 11)
add_edge!(batoidea, 9, 10)
add_edge!(batoidea, 10, 11)

filter_fn = make_unweighted_filter(cross, cross_boundary)
result = filter_fn(batoidea, nothing, cross_boundary)

println("BATOIDEA is a valid crossing gadget replacement: ", result !== nothing)
println("constant_offset = ", result.constant_offset)   # expected: 2.0

# We can also use is_gadget_replacement for direct pairwise checks:
valid, offset = is_gadget_replacement(cross, batoidea, cross_boundary, cross_boundary)
println("is_gadget_replacement: valid=$(valid), offset=$(offset)")

# ============================================================
# Step 3: Search a small triangular UDG dataset
# ============================================================
#
# Load the pre-generated 3×3 triangular UDG dataset (12 graphs, up to 9
# vertices each).  Boundary pins are always at vertex indices [1,2,3,4].
# The JIT is already warm from Step 2, so this completes quickly.

data_path = pkgdir(GadgetSearch, "data", "grid_udgs", "m3n3pad1_min3max9_direct4.g6")
loader = GraphLoader(data_path; pinset=cross_boundary)

results = search_unweighted_gadgets(cross, cross_boundary, loader)
println("Replacements found in m3n3 dataset: $(length(results))")

# ============================================================
# Step 4: Full search over all pre-generated datasets  [offline]
# ============================================================
#
# The block below is only included when running the script directly
# (it is stripped from the documentation build via the #src tag).
# Full results on the m4n4 dataset: 26 replacements found,
# vertices range 12-19, constant offsets 3.0 or 4.0.

#src datasets = [
#src     pkgdir(GadgetSearch, "data", "grid_udgs", "m3n3pad1_min3max9_direct4.g6"),
#src     pkgdir(GadgetSearch, "data", "grid_udgs", "m3n4pad1_min3max12_direct4.g6"),
#src     pkgdir(GadgetSearch, "data", "grid_udgs", "m3n5pad1_min3max15_direct4.g6"),
#src     pkgdir(GadgetSearch, "data", "grid_udgs", "m4n4pad1_min4max16_direct4.g6"),
#src ]
#src
#src all_results = UnweightedGadget[]
#src for path in datasets
#src     isfile(path) || (@warn "Dataset not found, skipping: $path"; continue)
#src     loader_i = GraphLoader(path; pinset=cross_boundary)
#src     @info "Searching $(basename(path))  ($(length(loader_i)) graphs)"
#src     t = @elapsed r = search_unweighted_gadgets(cross, cross_boundary, loader_i)
#src     @info "  → found $(length(r)) replacement(s) in $(round(t; digits=2))s"
#src     append!(all_results, r)
#src end
#src
#src println("Total crossing gadget replacements found: $(length(all_results))")
#src for (i, r) in enumerate(all_results)
#src     valid, offset = is_gadget_replacement(cross, r.replacement_graph,
#src                                            cross_boundary, r.boundary_vertices)
#src     println("Gadget #$(i): $(nv(r.replacement_graph))v/$(ne(r.replacement_graph))e  offset=$(r.constant_offset)  valid=$(valid)")
#src     out_path = joinpath(pkgdir(GadgetSearch, "examples"), "crossing_gadget_$(i).svg")
#src     GadgetSearch.plot_graph(r.replacement_graph, out_path; pos=r.pos,
#src                             plot_size=400, margin=40, vertex_size=12,
#src                             vertex_label_size=14, edge_width=2)
#src end
