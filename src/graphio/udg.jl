"""
    unit_disk_graph(locs::AbstractVector, unit::Real) -> SimpleGraph

Create a unit disk graph from given locations. Two vertices are connected 
if their Euclidean distance is less than the unit distance.

# Arguments
- `locs`: Vector of vertex locations
- `unit`: Unit distance threshold for edge creation

# Returns
- `SimpleGraph`: The constructed unit disk graph
"""
function unit_disk_graph(locs::AbstractVector, unit::Real)
    n = length(locs)
    g = SimpleGraph(n)
    for i = 1:n, j = i+1:n
        if sum(abs2, locs[i] .- locs[j]) < unit ^ 2
            add_edge!(g, i, j)
        end
    end
    return g
end

"""
    build_triangular_graph(grid_coords::Vector{Tuple{Int,Int}}) -> SimpleGraph{Int}

Build a `SimpleGraph` from integer grid coordinates using triangular-lattice
adjacency.  This is the shared building block for `generate_triangular_udg_subsets`,
`dedup_inner_subsets`, and downstream scripts.
"""
function build_triangular_graph(grid_coords::Vector{Tuple{Int,Int}})
    n = length(grid_coords)
    g = SimpleGraph(n)
    for a in 1:n, b in a+1:n
        ia, ja = grid_coords[a]
        ib, jb = grid_coords[b]
        triangular_adjacency(ia, ja, ib, jb) && add_edge!(g, a, b)
    end
    return g
end

"""
    generate_full_grid_graph(lattice::LatticeType, nx::Int, ny::Int; path::String="grid.g6") -> String

Generate complete graphs on a grid lattice (without boundary expansion) and save to file.
All vertices are connected to each other, regardless of distance.

# Arguments  
- `lattice`: Type of lattice (Square or Triangular) - determines physical positions
- `nx`: Number of grid points in x direction
- `ny`: Number of grid points in y direction
- `path`: Output file path for saving graphs (default: "grid.g6")

# Returns
- `String`: Path to the saved graph file

# Details
Generates a single complete graph on the nx×ny grid. The physical positions
follow the lattice geometry, but edges connect all pairs of vertices.
"""
function generate_full_grid_graph(lattice::LatticeType, nx::Int, ny::Int; path::String="grid.g6")
    # grid points (no boundary expansion)
    grid_points = get_physical_positions(lattice, 
        vec(Tuple{Int, Int}[(x, y) for x in 1:nx, y in 1:ny]))

    n = length(grid_points)
    g = Graphs.complete_graph(n)

    results = Tuple{SimpleGraph{Int}, Vector{Tuple{Float64, Float64}}}[(g, grid_points)]
    
    save_graph(results, path)
    @info "Generated complete graph with $n vertices on $(nx)×$(ny) grid"
    return path
end

"""
    generate_full_grid_udg(lattice::LatticeType, nx::Int, ny::Int; path::String="udg.g6") -> String

Generate unit disk graphs on a grid lattice with four boundary pins and save to file.

# Arguments  
- `lattice`: Type of lattice (Square or Triangular)
- `nx`: Number of inner positions in x direction
- `ny`: Number of inner positions in y direction
- `path`: Output file path for saving graphs (default: "udg.g6")

# Returns
- `String`: Path to the saved graph file

# Details
Generates all possible UDGs by placing pins on boundary positions and 
connecting vertices within unit distance on the specified lattice type.
"""
function generate_full_grid_udg(lattice::LatticeType, nx::Int, ny::Int; path::String="udg.g6")
    # get possible pin positions on the boundary
    top_candidates, bottom_candidates, left_candidates, right_candidates = get_pin_positions(lattice, nx, ny)

    # inner points are fixed
    inner_points = get_inner_points(lattice, nx, ny)

    radius = get_radius(lattice)

    results = Tuple{SimpleGraph{Int}, Vector{Tuple{Float64, Float64}}}[]

    for (top, bottom, left, right) in Iterators.product(
            top_candidates, bottom_candidates, left_candidates, right_candidates)

        selected = vcat([top, right, bottom, left], inner_points)

        g = unit_disk_graph(selected, radius)

        push!(results, (g, selected))
    end
    @info "pinset in generated graphs: [1,2,3,4]"
    return _process_and_save_graphs(results, path)
end

"""
    generate_triangular_udg_subsets(nx::Int, ny::Int; min_inner::Int=3, max_inner::Int=nx*ny, path::String="tri_subsets.g6") -> String

Generate triangular-lattice Unit Disk Graphs by enumerating **all subsets** of
inner grid positions, analogous to the pre-generated KSG (`grid_udgs`) datasets.

# Arguments
- `nx`: Number of inner columns
- `ny`: Number of inner rows
- `min_inner`: Minimum number of inner vertices to include (default 3)
- `max_inner`: Maximum number of inner vertices to include (default `nx*ny`)
- `path`: Output file path (default `"tri_subsets.g6"`)

# Details
The full padded grid has size `(nx+2) × (ny+2)`.  
Four **fixed** boundary pins are placed at the midpoints of the four sides:
- pin 1 (top):    `(⌈(nx+2)/2⌉, 1)`
- pin 2 (right):  `(nx+2, ⌈(ny+2)/2⌉)`
- pin 3 (bottom): `(⌈(nx+2)/2⌉, ny+2)`
- pin 4 (left):   `(1, ⌈(ny+2)/2⌉)`

Edges are determined by `triangular_adjacency` on integer `(i,j)` grid
coordinates (no floating-point distance computation). Physical positions use
the parallelogram layout `(x + (y−1)·0.5, y·√3/2)` for visualization only.

Vertex indices `[1,2,3,4]` in every generated graph correspond to the four
boundary pins, matching the `pinset=[1,2,3,4]` convention of `GraphLoader`.

# Returns
- `String`: Path to the saved graph file
"""
function generate_triangular_udg_subsets(nx::Int, ny::Int;
                                          min_inner::Int = 3,
                                          max_inner::Int = nx * ny,
                                          path::String = "tri_subsets.g6")
    Mx = nx + 2
    My = ny + 2

    # Fixed boundary pins at midpoints of the four sides
    cx = (Mx + 1) ÷ 2
    cy = (My + 1) ÷ 2

    pin_grid = Tuple{Int,Int}[(cx, 1), (Mx, cy), (cx, My), (1, cy)]
    pin_phys = get_physical_positions(Triangular(), pin_grid)

    # All inner grid positions
    inner_grid = Tuple{Int,Int}[(x, y) for x in 2:Mx-1 for y in 2:My-1]
    inner_phys = get_physical_positions(Triangular(), inner_grid)

    n_inner = length(inner_grid)
    max_k   = min(max_inner, n_inner)

    results = Tuple{SimpleGraph{Int}, Vector{Tuple{Float64,Float64}}}[]

    for k in min_inner:max_k
        for subset in combinations(1:n_inner, k)
            # Integer grid coords for all vertices: [pins..., selected inner...]
            all_grid = vcat(pin_grid, inner_grid[subset])
            all_phys = vcat(pin_phys, inner_phys[subset])
            n = length(all_grid)

            g = build_triangular_graph(all_grid)

            push!(results, (g, all_phys))
        end
    end

    @info "Generated $(length(results)) triangular subset UDGs ($(nx)×$(ny) inner grid, inner vertices $(min_inner)–$(max_k))"
    @info "pinset in generated graphs: [1,2,3,4]"
    return _process_and_save_graphs(results, path)
end

"""
    dedup_inner_subsets(inner_grid::Vector{Tuple{Int,Int}}, k::Int) -> Vector{Vector{Int}}

Enumerate all `k`-subsets of `inner_grid`, build the corresponding inner-only
triangular-lattice subgraphs, and keep one representative per graph-isomorphism
class.

If `shortg` is available, it is used for exact canonical deduplication.
Otherwise, the function falls back to deduplicating by graph6 strings, which is
still effective but may miss some isomorphic classes when vertex order differs.
"""
function dedup_inner_subsets(inner_grid::Vector{Tuple{Int,Int}}, k::Int)
    n = length(inner_grid)
    k > n && return Vector{Int}[]

    all_subsets = collect(Combinatorics.combinations(1:n, k))
    isempty(all_subsets) && return Vector{Int}[]

    graphs = [build_triangular_graph(inner_grid[subset]) for subset in all_subsets]

    if _has_shortg()
        temp_file = tempname()
        mapping_file = tempname()
        try
            open(temp_file, "w") do io
                for g in graphs
                    println(io, GraphIO.Graph6._graphToG6String(g)[11:end])
                end
            end

            _call_shortg(temp_file, mapping_file)
            canon2orig, _ = _parse_shortg_mapping(mapping_file)
            rep_indices = sort!([first(v) for v in values(canon2orig)])
            return [all_subsets[i] for i in rep_indices]
        finally
            isfile(temp_file) && rm(temp_file; force=true)
            isfile(mapping_file) && rm(mapping_file; force=true)
        end
    end

    @warn "shortg not available; falling back to g6-string dedup (less thorough)"
    seen = Dict{String, Int}()
    rep_indices = Int[]
    for (i, g) in enumerate(graphs)
        g6 = GraphIO.Graph6._graphToG6String(g)
        if !haskey(seen, g6)
            seen[g6] = i
            push!(rep_indices, i)
        end
    end
    return [all_subsets[i] for i in rep_indices]
end

