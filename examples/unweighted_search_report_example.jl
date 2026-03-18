using GadgetSearch
using Graphs

function cross_graph()
    g = SimpleGraph(4)
    add_edge!(g, 1, 3)
    add_edge!(g, 2, 4)
    return g
end

function subdivide_edge(g::SimpleGraph, u::Int, v::Int)
    h = copy(g)
    rem_edge!(h, u, v)
    add_vertex!(h)
    w = nv(h)
    add_edge!(h, u, w)
    add_edge!(h, w, v)
    return h
end

target_graph = cross_graph()
target_boundary = [1, 2, 3, 4]
sample_candidate = subdivide_edge(subdivide_edge(target_graph, 1, 3), 2, 4)

outpath = pkgdir(GadgetSearch, "examples", "unweighted_search_report.svg")

plot_unweighted_search_report(
    target_graph,
    target_boundary,
    outpath;
    max_added_vertices=1,
    preserve_boundary_roles=true,
    sample_candidate=sample_candidate,
    sample_candidate_boundary=target_boundary,
    prefilter=true,
)

println("Saved report visualization to $outpath")
