using GadgetSearch

const EXAMPLES_DIR = pkgdir(GadgetSearch, "examples")

triangular_positions(nx::Int, ny::Int) =
    GadgetSearch.get_physical_positions(Triangular(), vec(Tuple{Int, Int}[(i, j) for i in 1:nx, j in 1:ny]))

function plot_lattice_examples()
    full_graph = triangular_lattice_graph(4, 4)
    full_positions = triangular_positions(4, 4)
    full_path = joinpath(EXAMPLES_DIR, "triangular_lattice_4x4.svg")
    GadgetSearch.plot_graph(full_graph, full_path; pos=full_positions, plot_size=700, vertex_size=8, vertex_label_size=12)

    small_graph = triangular_lattice_graph(3, 3)
    small_positions = triangular_positions(3, 3)
    small_path = joinpath(EXAMPLES_DIR, "triangular_lattice_3x3.svg")
    GadgetSearch.plot_graph(small_graph, small_path; pos=small_positions, plot_size=600, vertex_size=10, vertex_label_size=14)

    return (full_path, small_path)
end

function ensure_subset_dataset()
    dataset_path = joinpath(EXAMPLES_DIR, "triangular_subset_dataset.g6")
    if !isfile(dataset_path)
        generate_triangular_udg_subsets(3, 3; subset_sizes=2:3, deduplicate=true, path=dataset_path)
    end
    return dataset_path
end

function plot_subset_examples(dataset_path::String; nplots::Int=3)
    loader = GraphLoader(dataset_path)
    plot_count = min(nplots, length(loader))
    saved_paths = String[]

    for idx in 1:plot_count
        graph = loader[idx]
        positions = loader.layout[idx]
        positions === nothing && continue

        outpath = joinpath(EXAMPLES_DIR, "triangular_subset_$(idx).svg")
        GadgetSearch.plot_graph(graph, outpath; pos=positions, plot_size=500, vertex_size=12, vertex_label_size=16)
        push!(saved_paths, outpath)
    end

    return saved_paths
end

function main()
    lattice_paths = plot_lattice_examples()
    dataset_path = ensure_subset_dataset()
    subset_paths = plot_subset_examples(dataset_path)

    println("Triangular lattice plots:")
    foreach(println, lattice_paths)
    println("Subset dataset:")
    println(dataset_path)
    println("Subset plots:")
    foreach(println, subset_paths)
end

main()

