# # Search for Rydberg Gadgets on Triangular Lattice
#
# This example demonstrates how to search for **Rydberg gadgets** using:
# - **RydbergModel**: State space is restricted to Maximal Independent Sets (MIS)
# - **Unit Disk Graphs (UDG)**: Edges connect vertices within unit distance
# - **Vertex weights only**: Energy E(σ) = Σᵢ hᵢσᵢ
#
# Notes:
# - The search can be expensive. Tweak `max_samples`, lattice size, and `pin_candidates`

using GadgetSearch
using HiGHS
using Combinatorics
using FileIO, ImageShow

# Define truth tables for 2-input, 1-output logic gates
truth_tables = [
    TruthTableConstraint(BitMatrix([0 0 0; 1 0 1; 0 1 1; 1 1 1])),   # OR
    TruthTableConstraint(BitMatrix([0 0 0; 1 0 0; 0 1 0; 1 1 1])),   # AND
    TruthTableConstraint(BitMatrix([0 0 1; 1 0 1; 0 1 1; 1 1 0])),   # NAND
    TruthTableConstraint(BitMatrix([0 0 1; 1 0 0; 0 1 0; 1 1 0])),   # NOR
    TruthTableConstraint(BitMatrix([0 0 0; 1 0 1; 0 1 1; 1 1 0]))    # XOR
]

# Generate Unit Disk Graph dataset on triangular lattice
# The UDG constraint means only nearby atoms can interact (Rydberg blockade)
generate_full_grid_udg(Triangular(), 2, 2; path=pkgdir(GadgetSearch, "examples", "rydberg_dataset.g6"))

dataloader = GraphLoader(pkgdir(GadgetSearch, "examples", "rydberg_dataset.g6"))

# Search using RydbergModel explicitly
# - RydbergModel uses MIS (Maximal Independent Sets) as the state space
# - Only vertex weights (h) are optimized
# - Objective minimizes sum of vertex weights
results, failed = search_gadgets(
    RydbergModel,           # Explicitly use Rydberg model
    dataloader, 
    truth_tables;
    optimizer=HiGHS.Optimizer, 
    pin_candidates=collect(Combinatorics.combinations(1:4, 3)), 
    allow_defect=true,  
    objective=h -> sum(h),  # Only vertex weights for Rydberg
    save_path=joinpath(pkgdir(GadgetSearch, "examples"), "triangular_Rydberg_results.json"), 
    max_result_num=10,
    max_samples=10000,
    check_connectivity=true
)

println("================================")
println("   RYDBERG MODEL SEARCH RESULTS")
println("================================")
@info "Results: $(length(results)) truth tables searched"
for (i, tt_results) in enumerate(results)
    @info "Truth table $(i-1): $(length(tt_results)) gadgets found"
end

@info "Cache statistics: $(GadgetSearch.get_cache_stats())"
GadgetSearch.clear_cache!()

# ## Visualize found gadgets

labels = ["OR", "AND", "NAND", "NOR", "XOR"]

for (i, label) in enumerate(labels)
    if !isempty(results[i])
        outpath = pkgdir(GadgetSearch, "examples", "gadget_$(label).png")
        GadgetSearch.plot_gadget(results[i][1], outpath; show_weights=true, round_weights=true)
        @info "Saved $label gate to $outpath"
    end
end

# ## Check correctness using Rydberg model
# The check_gadget function validates ground states using MIS

for (i, label) in enumerate(labels)
    if !isempty(results[i])
        gadget = results[i][1]
        @info """===== Checking $(label) gate (Rydberg/MIS) =====
        Pins: $(gadget.pins)
        Vertex weights: $(gadget.vertex_weights)
        $(check_gadget_rydberg(gadget; _return_info=true))
        """
    end
end
