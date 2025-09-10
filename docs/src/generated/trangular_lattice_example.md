```@meta
EditURL = "../../../examples/trangular_lattice_example.jl"
```

# Search for gadgets on triangular lattice

````@example trangular_lattice_example
using Pkg
Pkg.add(["Combinatorics", "HiGHS"])

using GadgetSearch
using HiGHS
using Combinatorics
````

truth_table = GadgetSearch.generic_rule(110, (3, 1))

````@example trangular_lattice_example
truth_table = BitMatrix.([
    [0 0 0; 1 0 1; 0 1 1; 1 1 1],   # OR
    [0 0 0; 1 0 0; 0 1 0; 1 1 1],   # AND
    [0 0 1; 1 0 1; 0 1 1; 1 1 0],   # NAND = not(AND)
    [0 0 1; 1 0 0; 0 1 0; 1 1 0],   # NOR  = not(OR)
    [0 0 0; 1 0 1; 0 1 1; 1 1 0]    # XOR
]);

generate_full_grid_udg(Triangular(), 2, 2; path="dataset.g6")

dataloader = GraphLoader("dataset.g6")
````

Performance optimization: Reduce max_samples for faster initial search
Increase max_samples (e.g., 1000) for thorough search if needed

````@example trangular_lattice_example
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
       )

@show size(results)
````

Display results structure

````@example trangular_lattice_example
println("--------------------------------")
@info "Results structure: $(length(results)) truth tables"
for (i, tt_results) in enumerate(results)
    @info "Truth table $(i-1): $(length(tt_results)) gadgets found"
end
````

Check cache performance

````@example trangular_lattice_example
@info "Cache statistics after search: $(GadgetSearch.get_cache_stats())"
````

Uncomment to clear cache if memory is an issue:
GadgetSearch.clear_cache!()

Plot gadgets for each truth table

````@example trangular_lattice_example
for (tt_idx, tt_results) in enumerate(results)
    for (gadget_idx, gadget) in enumerate(tt_results)
        GadgetSearch.plot_gadget(gadget, "gadget$(tt_idx)_$(gadget_idx).pdf"; show_weights=true)
    end
end
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

