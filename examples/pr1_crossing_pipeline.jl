# # PR1 Crossing Pipeline (line-by-line runnable)
#
# This script demonstrates three parts of the unweighted crossing workflow:
# 1) equivalent_representations
# 2) logical flip utilities
# 3) search_unweighted_gadgets
#
# It is intentionally organized into small, single-purpose functions so users can
# execute each section line by line in the REPL and inspect output immediately.

using GadgetSearch
using Graphs
using GenericTensorNetworks: content

const OUTPUT_DIR = pkgdir(GadgetSearch, "examples", "pr1_pipeline_output")

function build_cross_graph()
    g = SimpleGraph(4)
    add_edge!(g, 1, 3)
    add_edge!(g, 2, 4)
    return g
end

function subdivide_edge(g::SimpleGraph{Int}, u::Int, v::Int, count::Int=1)
    count >= 0 || throw(ArgumentError("count must be non-negative"))
    has_edge(g, u, v) || throw(ArgumentError("edge ($u, $v) does not exist"))
    count == 0 && return copy(g)

    expanded = copy(g)
    rem_edge!(expanded, u, v)
    previous = u
    for _ in 1:count
        add_vertex!(expanded)
        next_vertex = nv(expanded)
        add_edge!(expanded, previous, next_vertex)
        previous = next_vertex
    end
    add_edge!(expanded, previous, v)
    return expanded
end

function subdivide_cross_graph(specs::Vector{Tuple{Int, Int, Int}})
    g = build_cross_graph()
    for (u, v, count) in specs
        g = subdivide_edge(g, u, v, count)
    end
    return g
end

function print_graph_summary(name::AbstractString, g::SimpleGraph{Int})
    println("[$name] vertices=$(nv(g)), edges=$(ne(g)), edge_list=$(collect(edges(g)))")
end

function ensure_output_dir()
    mkpath(OUTPUT_DIR)
    println("Output directory: $OUTPUT_DIR")
end

function plot_canonical_crossing(g::SimpleGraph{Int})
    path = joinpath(OUTPUT_DIR, "crossing_canonical.svg")
    GadgetSearch.plot_graph(g, path)
    println("Saved canonical crossing plot -> $path")
    return path
end

function plot_single_subdivision(g::SimpleGraph{Int})
    path = joinpath(OUTPUT_DIR, "crossing_subdivision_single.svg")
    GadgetSearch.plot_graph(g, path)
    println("Saved single-subdivision plot -> $path")
    return path
end

function plot_double_subdivision(g::SimpleGraph{Int})
    path = joinpath(OUTPUT_DIR, "crossing_subdivision_double.svg")
    GadgetSearch.plot_graph(g, path)
    println("Saved double-subdivision plot -> $path")
    return path
end

function plot_found_replacement(g::SimpleGraph{Int})
    path = joinpath(OUTPUT_DIR, "crossing_search_match.svg")
    GadgetSearch.plot_graph(g, path)
    println("Saved search-match plot -> $path")
    return path
end

function plot_gadget_batch(canonical::SimpleGraph{Int}, once_subdivided::SimpleGraph{Int}, twice_subdivided::SimpleGraph{Int})
    println("\n=== Batch plotting ===")
    return (
        canonical=plot_canonical_crossing(canonical),
        once=plot_single_subdivision(once_subdivided),
        twice=plot_double_subdivision(twice_subdivided),
    )
end

function run_equivalent_representations_demo(target_graph::SimpleGraph{Int}, target_boundary::Vector{Int})
    println("\n=== Module: equivalent_representations ===")
    reprs = equivalent_representations(target_graph, target_boundary; max_added_vertices=2)
    println("Representations found: $(length(reprs))")
    for (i, (g, boundary)) in enumerate(reprs)
        println("  repr[$i]: vertices=$(nv(g)), edges=$(ne(g)), boundary=$boundary")
    end
    return reprs
end

function run_flip_demo(target_graph::SimpleGraph{Int}, target_boundary::Vector{Int})
    println("\n=== Module: flip ===")
    reduced = Float64.(content.(calculate_reduced_alpha_tensor(target_graph, target_boundary)))
    patterns = generate_flip_patterns(length(target_boundary))
    println("Flip patterns: $(length(patterns))")

    for (mask, desc) in patterns
        flipped = apply_flip_to_tensor(reduced, mask)
        finite_deltas = [f - b for (f, b) in zip(vec(flipped), vec(reduced)) if isfinite(f) && isfinite(b)]
        example_delta = ""
        if isempty(finite_deltas)
            example_delta = "n/a"
        else
            example_delta = string(first(finite_deltas))
        end
        println("  $desc, mask=$mask, finite_delta_example=$example_delta")
    end
end

function build_crossing_demo_candidates(target_graph::SimpleGraph{Int})
    candidate_a = target_graph
    candidate_b = subdivide_cross_graph([(1, 3, 1)])
    candidate_c = subdivide_cross_graph([(1, 3, 1), (2, 4, 1)])
    return [candidate_a, candidate_b, candidate_c]
end

function build_loader_from_candidates(candidate_graphs::Vector{SimpleGraph{Int}}, target_boundary::Vector{Int})
    isempty(candidate_graphs) && throw(ArgumentError("candidate_graphs must be non-empty"))
    dataset = GraphDataset(graph_to_g6.(candidate_graphs))
    return GraphLoader(dataset, pinset=target_boundary)
end

function run_search_demo(target_graph::SimpleGraph{Int}, target_boundary::Vector{Int}, candidate_graphs::Vector{SimpleGraph{Int}})
    println("\n=== Module: search ===")
    loader = build_loader_from_candidates(candidate_graphs, target_boundary)

    results = search_unweighted_gadgets(
        target_graph,
        target_boundary,
        loader;
        max_added_vertices=2,
        include_logical_flips=true,
        max_results=10,
    )

    println("Search hits: $(length(results))")
    for (i, result) in enumerate(results)
        println("  hit[$i]: boundary=$(result.boundary_vertices), offset=$(result.constant_offset), vertices=$(nv(result.replacement_graph))")
    end
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("PR1 crossing pipeline demo (line-by-line friendly)")
    ensure_output_dir()

    target_graph = build_cross_graph()
    target_boundary = [1, 2, 3, 4]

    println("\n=== Base graphs ===")
    once_subdivided = subdivide_cross_graph([(1, 3, 1)])
    twice_subdivided = subdivide_cross_graph([(1, 3, 1), (2, 4, 1)])
    print_graph_summary("canonical", target_graph)
    print_graph_summary("subdivided_once", once_subdivided)
    print_graph_summary("subdivided_twice", twice_subdivided)

    plot_gadget_batch(target_graph, once_subdivided, twice_subdivided)
    reprs = run_equivalent_representations_demo(target_graph, target_boundary)
    run_flip_demo(target_graph, target_boundary)
    candidates = build_crossing_demo_candidates(target_graph)
    results = run_search_demo(target_graph, target_boundary, candidates)

    if !isempty(results)
        plot_found_replacement(results[1].replacement_graph)
    end

    println("\nDone. You can now inspect outputs in: $OUTPUT_DIR")
    println("Tip: execute function calls above one section at a time in REPL.")
end
