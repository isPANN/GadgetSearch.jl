# Search for Crossing Gadget (QUBO) on Triangular Lattice
#
# pin1 → pin3 (wire 1), pin2 → pin4 (wire 2)
# opposite_pairs = [(1,3), (2,4)]: each wire's pins on opposite boundary sides

using GadgetSearch
using Gurobi

crossing_constraint = TruthTableConstraint(Bool[
    0 0 0 0
    0 1 0 1
    1 0 1 0
    1 1 1 1
])

opposite_pairs = [(1, 3), (2, 4)]

# Generate UDG candidates on triangular lattice
nx, ny = 3, 3
dataset_path = joinpath(pkgdir(GadgetSearch), "examples", "crossing_triangular_dataset.jsonl")
generate_udg_with_open_pins(Triangular(), nx, ny; pin_count=4, path=dataset_path)

loader = GraphLoader(dataset_path)
println("$(length(loader)) candidate graphs")

env = Gurobi.Env()

results, failed = search_gadgets(
    QUBOModel, loader, [crossing_constraint];
    optimizer=Gurobi.Optimizer, env=env,
    allow_defect=true,
    max_result_num=5,
    max_samples=5000,
    check_connectivity=true,
    opposite_pairs=opposite_pairs,
    objective=(h, J) -> sum(h .^ 2) + sum(J .^ 2),
    save_path=joinpath(pkgdir(GadgetSearch), "examples", "crossing_triangular_results.json")
)

@show failed

if !isempty(results[1])
    for (i, gadget) in enumerate(results[1])
        println("\n===== Crossing Gadget #$i =====")
        println("Pins: ", gadget.pins)
        println("Vertex weights: ", gadget.vertex_weights)
        println("Edge weights: ", gadget.edge_weights)
        println("Edge list: ", gadget.edge_list)
        outpath = joinpath(pkgdir(GadgetSearch), "examples", "crossing_triangular_gadget_$i.png")
        plot_gadget(gadget, outpath;
            show_weights=true, show_edge_weights=true, round_weights=true)
    end
    println("\n===== Validation =====")
    println(check_gadget_qubo(results[1][1]; _return_info=true))
else
    println("No crossing gadget found. Try increasing nx, ny.")
end

GadgetSearch.clear_cache!()
