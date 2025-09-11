# # Search for gadgets on triangular lattice
#
# This example demonstrates how to:
# - Search for gadgets logically equivalent to several 2-input, 1-output logic gates (OR/AND/NAND/NOR/XOR) on the triangular-lattice
# - Visualize found gadgets
# Please check [our paper (to be published)] for the searching details.
#
# Notes:
# - The search can be expensive. Tweak the choices of `max_samples`, lattice size, and `pin_candidates`

using GadgetSearch
using HiGHS
using Combinatorics
using FileIO, ImageShow

# We can simply specify the truth tables of the logic gates we want to search for.

truth_table = BitMatrix.([
    [0 0 0; 1 0 1; 0 1 1; 1 1 1],   # OR
    [0 0 0; 1 0 0; 0 1 0; 1 1 1],   # AND
    [0 0 1; 1 0 1; 0 1 1; 1 1 0],   # NAND = not(AND)
    [0 0 1; 1 0 0; 0 1 0; 1 1 0],   # NOR  = not(OR)
    [0 0 0; 1 0 1; 0 1 1; 1 1 0]    # XOR
]);

# Before we start the search, we need to generate the candidate UDG dataset.
# Here we generate a 2x2 triangular-lattice UDG dataset. Here we add four boundary pins to the lattice.

generate_full_grid_udg(Triangular(), 2, 2; path=pkgdir(GadgetSearch, "examples", "dataset.g6"))

# Then we load the dataset.

dataloader = GraphLoader(pkgdir(GadgetSearch, "examples", "dataset.g6"))

# Begin the search process:
# - Because we add four boundary pins to the lattice, we must specify both which pins are used and their order in `pin_candidates`.
# - If `allow_defect` is true, a vertex may be assigned weight 0, meaning the vertex and its associated edges can effectively be removed from the gadget.
# - The objective function minimizes the sum of vertex weights, encouraging the simplest possible gadget.
# - `max_result_num` limits how many gadgets are found for each truth table; once this limit is reached, the search for that truth table stops.
# - For each graph under a given pin assignment, the ground-state configurations may still have many possibilities, leading to a large search space. We use `max_samples` to enable early termination, though this typically results in incomplete exploration.
# - `check_connectivity` ensures that the resulting gadget remains connected.

results, failed = search_by_truth_tables(
           dataloader, 
           truth_table;
           optimizer=HiGHS.Optimizer, 
           pin_candidates=collect(Combinatorics.combinations(1:4, 3)), 
           allow_defect=true,  
           objective=x->sum(x), 
           max_result_num=10,
           max_samples=10000,
           check_connectivity=true
       );

println("--------------------------------")
@info "Results structure: $(length(results)) truth tables"
for (i, tt_results) in enumerate(results)
    @info "Truth table $(i-1): $(length(tt_results)) gadgets found"
end

# Check cache performance

@info "Cache statistics after search: $(GadgetSearch.get_cache_stats())"

# Clear cache if memory is an issue:

GadgetSearch.clear_cache!()

# ## Gadget examples on triangular lattice

# Note the gadget for the same truth table is not unique. 

# ### OR gate

outpath = pkgdir(GadgetSearch, "examples", "gadget_OR.png")
GadgetSearch.plot_gadget(results[1][1], outpath; show_weights=true, round_weights=true)
load(outpath)

# ### AND gate

outpath = pkgdir(GadgetSearch, "examples", "gadget_AND.png")
GadgetSearch.plot_gadget(results[2][1], outpath; show_weights=true, round_weights=true)
load(outpath)

# ### NAND gate

outpath = pkgdir(GadgetSearch, "examples", "gadget_NAND.png")
GadgetSearch.plot_gadget(results[3][1], outpath; show_weights=true, round_weights=true)
load(outpath)

# ### NOR gate

outpath = pkgdir(GadgetSearch, "examples", "gadget_NOR.png")
GadgetSearch.plot_gadget(results[4][1], outpath; show_weights=true, round_weights=true)
load(outpath)

# ### XOR gate

outpath = pkgdir(GadgetSearch, "examples", "gadget_XOR.png")
GadgetSearch.plot_gadget(results[5][1], outpath; show_weights=true, round_weights=true)
load(outpath)

# ### Check correctness
# We can use `check_gadget` to check the correctness, i.e. the ground states, of the found gadgets.

labels = ["OR", "AND", "NAND", "NOR", "XOR"]
for i in 1:5
    @info """===== Checking correctness of $(labels[i]) =====
    => Truth table:
    $(truth_table[i])
    => Found gadget:
    $(GadgetSearch.check_gadget(results[i][1]; _return_info=true))
    """
end