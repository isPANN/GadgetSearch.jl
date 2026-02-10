# # Search for Unweighted Rydberg Gadgets on Triangular Lattice
#
# This example demonstrates how to search for **unweighted Rydberg gadgets** using:
# - **RydbergUnweightedModel**: All vertex weights are fixed to 1 (uniform)
# - **Unit Disk Graphs (UDG)**: Edges connect vertices within unit distance
# - **No optimizer needed**: Gadgets are found purely by graph structure
#
# Reference: Liu et al., "Computer-assisted gadget design and problem reduction
# of unweighted maximum independent set"
#
# In the unweighted model, ground states are the Maximum Independent Sets
# (MIS with the largest cardinality). The search checks if the pin patterns
# of these maximum MIS naturally match the desired truth table.

using GadgetSearch
using Combinatorics

# Define truth tables for logic gates (same as weighted Rydberg example)
truth_tables = [
    TruthTableConstraint(BitMatrix([0 0 0; 1 0 1; 0 1 1; 1 1 1])),   # OR
    TruthTableConstraint(BitMatrix([0 0 0; 1 0 0; 0 1 0; 1 1 1])),   # AND
    TruthTableConstraint(BitMatrix([0 0 1; 1 0 1; 0 1 1; 1 1 0])),   # NAND
    TruthTableConstraint(BitMatrix([0 0 1; 1 0 0; 0 1 0; 1 1 0])),   # NOR
    TruthTableConstraint(BitMatrix([0 0 0; 1 0 1; 0 1 1; 1 1 0]))    # XOR
]

# Generate Unit Disk Graph dataset on triangular lattice
generate_full_grid_udg(Triangular(), 3, 3; path=pkgdir(GadgetSearch, "examples", "unweighted_dataset.g6"))

dataloader = GraphLoader(pkgdir(GadgetSearch, "examples", "unweighted_dataset.g6"))

# Search using RydbergUnweightedModel
# Key differences from RydbergModel:
# - No optimizer needed (weights are all 1)
# - No objective function needed
# - Search is faster (pure feasibility check, no optimization)
results, failed = search_gadgets(
    RydbergUnweightedModel,  # Unweighted Rydberg model
    dataloader, 
    truth_tables;
    # No optimizer needed!
    pin_candidates=collect(Combinatorics.combinations(1:4, 3)), 
    save_path=joinpath(pkgdir(GadgetSearch, "examples"), "triangular_unweighted_results.json"), 
    max_result_num=5,
    max_samples=10000,
    check_connectivity=true
)

println("==========================================")
println("   RYDBERG UNWEIGHTED MODEL SEARCH RESULTS")
println("==========================================")
@info "Results: $(length(results)) truth tables searched"

labels = ["OR", "AND", "NAND", "NOR", "XOR"]
for (i, label) in enumerate(labels)
    @info "Truth table '$label': $(length(results[i])) gadgets found"
end

if !isempty(failed)
    @info "Failed constraints: $(length(failed))"
    @info "Note: Unweighted gadgets are much harder to find than weighted ones."
    @info "Try larger lattice sizes (e.g., 4x4 or 5x5) for more candidates."
end

@info "Cache statistics: $(GadgetSearch.get_cache_stats())"
GadgetSearch.clear_cache!()

# ## Visualize found gadgets (if any)

for (i, label) in enumerate(labels)
    if !isempty(results[i])
        outpath = pkgdir(GadgetSearch, "examples", "unweighted_gadget_$(label).png")
        GadgetSearch.plot_gadget(results[i][1], outpath; show_weights=true, round_weights=true)
        @info "Saved $label gate to $outpath"
    end
end

# ## Check correctness using unweighted model

for (i, label) in enumerate(labels)
    if !isempty(results[i])
        gadget = results[i][1]
        println("\n===== Checking $(label) gate (Rydberg Unweighted/MIS) =====")
        println("Pins: $(gadget.pins)")
        println("All vertex weights are 1: $(all(w -> w == 1.0, gadget.vertex_weights))")
        println(check_gadget_unweighted(gadget; _return_info=true))
    end
end

