using GadgetSearch
using Graphs

function cross_graph()
    g = SimpleGraph(4)
    add_edge!(g, 1, 3)
    add_edge!(g, 2, 4)
    return g
end

graph = cross_graph()
boundary = [1, 2, 3, 4]
outpath = pkgdir(GadgetSearch, "examples", "equivalent_representations_cross.svg")

reprs = equivalent_representations(graph, boundary; max_added_vertices=1)
println("Found $(length(reprs)) distinct equivalent representations.")

plot_equivalent_representations(
    graph,
    boundary,
    outpath;
    max_added_vertices=1,
    columns=3,
    panel_graph_size=150,
    show_reduced_tensor=false,
)

println("Saved visualization to $outpath")

