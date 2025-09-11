# # Search for gadgets on triangular lattice
using Pkg
Pkg.activate(pkgdir(GadgetSearch)*"/docs")

using GadgetSearch
using HiGHS
using Combinatorics
using FileIO

# truth_table = GadgetSearch.generic_rule(110, (3, 1))
truth_table = BitMatrix.([
    [0 0 0; 1 0 1; 0 1 1; 1 1 1],   # OR
    [0 0 0; 1 0 0; 0 1 0; 1 1 1],   # AND
    [0 0 1; 1 0 1; 0 1 1; 1 1 0],   # NAND = not(AND)
    [0 0 1; 1 0 0; 0 1 0; 1 1 0],   # NOR  = not(OR)
    [0 0 0; 1 0 1; 0 1 1; 1 1 0]    # XOR
]);

generate_full_grid_udg(Triangular(), 2, 2; path="dataset.g6")

dataloader = GraphLoader("dataset.g6")

# Performance optimization: Reduce max_samples for faster initial search
# Increase max_samples (e.g., 1000) for thorough search if needed
results, failed = search_by_truth_tables(
           dataloader, 
           truth_table;
           optimizer=HiGHS.Optimizer, 
           pin_candidates=collect(Combinatorics.combinations(1:4, 3)), 
           allow_defect=true,  
           objective=x->sum(x), 
           max_result_num=10,
           max_samples=10000,  # Reduced from default 100 for faster execution
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

# Plot one of the results for OR gate
GadgetSearch.plot_gadget(results[1][1], "gadget_OR.png"; show_weights=true, round_weights=true)
display("image/png", read("gadget_OR.png"))

# Plot one of the results for AND gate
GadgetSearch.plot_gadget(results[2][1], "gadget_AND.png"; show_weights=true, round_weights=true)
display("image/png", read("gadget_AND.png"))

# Plot one of the results for NAND gate
GadgetSearch.plot_gadget(results[3][1], "gadget_NAND.png"; show_weights=true, round_weights=true)
display("image/png", read("gadget_NAND.png"))

# Plot one of the results for NOR gate
GadgetSearch.plot_gadget(results[4][1], "gadget_NOR.png"; show_weights=true, round_weights=true)
display("image/png", read("gadget_NOR.png"))

# Plot one of the results for XOR gate
GadgetSearch.plot_gadget(results[5][1], "gadget_XOR.png"; show_weights=true, round_weights=true)
display("image/png", read("gadget_XOR.png"))

