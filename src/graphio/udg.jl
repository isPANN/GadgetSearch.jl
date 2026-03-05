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

abstract type LatticeType end
struct Square <: LatticeType end
struct Triangular <: LatticeType end

get_radius(::Square) = 1.5
get_radius(::Triangular) = 1.1

get_physical_positions(::Square, pos::Vector{Tuple{Int, Int}}) = Vector{Tuple{Float64, Float64}}(pos)
function get_physical_positions(::Triangular, pos::Vector{Tuple{Int, Int}})
    h = sqrt(3) / 2
    return Vector{Tuple{Float64, Float64}}([(x + (y - 1) * 0.5, y * h) for (x, y) in pos])
end

function get_pin_positions(lattice::LatticeType, nx::Int, ny::Int)
    top_candidates    = get_physical_positions(lattice, Tuple{Int, Int}[(1, y)     for y in 2:ny+1])
    bottom_candidates = get_physical_positions(lattice, Tuple{Int, Int}[(nx+2, y)  for y in 2:ny+1])
    left_candidates   = get_physical_positions(lattice, Tuple{Int, Int}[(x, 1)     for x in 2:nx+1])
    right_candidates  = get_physical_positions(lattice, Tuple{Int, Int}[(x, ny+2)  for x in 2:nx+1])
    return top_candidates, bottom_candidates, left_candidates, right_candidates
end

function get_inner_points(lattice::LatticeType, nx::Int, ny::Int)
    original_points = Tuple{Int, Int}[(x, y) for x in 2:nx+1, y in 2:ny+1]
    return get_physical_positions(lattice, vec(original_points))
end

"""
    triangular_adjacency(i1, j1, i2, j2) -> Bool

Return `true` if integer grid positions `(i1, j1)` and `(i2, j2)` are
nearest neighbours on the triangular lattice.

The triangular lattice is a rectangular grid with **one diagonal** added per
cell (the `(+1,−1)` / `(−1,+1)` direction). Under the parallelogram physical
mapping `(x, y) → (x + (y−1)·0.5,  y·√3/2)` all six resulting neighbours lie
at Euclidean distance exactly 1.0.

```
neighbours of (i,j):
  (i±1, j)   – horizontal
  (i, j±1)   – vertical
  (i+1, j−1) – diagonal ↘
  (i−1, j+1) – diagonal ↖
```

No floating-point computation is required.
"""
function triangular_adjacency(i1::Int, j1::Int, i2::Int, j2::Int)
    di, dj = i2 - i1, j2 - j1
    return (abs(di) == 1 && dj == 0)  ||  # horizontal
           (di == 0 && abs(dj) == 1)  ||  # vertical
           (di ==  1 && dj == -1)     ||  # diagonal ↘
           (di == -1 && dj ==  1)         # diagonal ↖
end

"""
    triangular_lattice_graph(nx, ny) -> (SimpleGraph{Int}, Vector{Tuple{Float64,Float64}})

Build the full `nx × ny` triangular-lattice graph directly from integer `(i, j)`
indices, **without** floating-point distance computation.

The connectivity rule is `triangular_adjacency`: each vertex connects to its
horizontal, vertical, and one diagonal neighbour (six neighbours in the interior).
Physical positions follow the parallelogram layout used by
`get_physical_positions(Triangular(), ...)` and are returned as the second
element of the tuple for use in visualisation or UDG construction.

# Arguments
- `nx`: Number of columns (x direction)
- `ny`: Number of rows    (y direction)

# Returns
- `(SimpleGraph{Int}, Vector{Tuple{Float64,Float64}})`: graph and vertex positions

# Example
```julia
g, pos = triangular_lattice_graph(3, 3)   # 9-vertex triangular lattice
```
"""
function triangular_lattice_graph(nx::Int, ny::Int)
    coords = Tuple{Int,Int}[(i, j) for i in 1:nx for j in 1:ny]
    n      = length(coords)
    g      = SimpleGraph(n)
    idx    = Dict{Tuple{Int,Int}, Int}(coords[k] => k for k in 1:n)

    for k in 1:n
        i, j = coords[k]
        for (di, dj) in ((1, 0), (0, 1), (1, -1))   # 3 undirected half-edges
            nb = (i + di, j + dj)
            if haskey(idx, nb)
                add_edge!(g, k, idx[nb])
            end
        end
    end

    pos = get_physical_positions(Triangular(), coords)
    return g, pos
end

"""
    complete_graph(n::Int) -> SimpleGraph

Create a complete graph with n vertices where every pair of vertices is connected.

# Arguments
- `n`: Number of vertices

# Returns
- `SimpleGraph`: The complete graph K_n
"""
function complete_graph(n::Int)
    g = SimpleGraph(n)
    for i = 1:n, j = i+1:n
        add_edge!(g, i, j)
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
    g = complete_graph(n)

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

    for top in top_candidates, bottom in bottom_candidates,
        left in left_candidates, right in right_candidates

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
    Ny_ = ny + 2

    # Fixed boundary pins at midpoints of the four sides
    cx = (Mx + 1) ÷ 2
    cy = (Ny_ + 1) ÷ 2

    pin_grid = Tuple{Int,Int}[(cx, 1), (Mx, cy), (cx, Ny_), (1, cy)]
    pin_phys = get_physical_positions(Triangular(), pin_grid)

    # All inner grid positions
    inner_grid = Tuple{Int,Int}[(x, y) for x in 2:Mx-1 for y in 2:Ny_-1]
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

            g = SimpleGraph(n)
            for a in 1:n, b in a+1:n
                ia, ja = all_grid[a]
                ib, jb = all_grid[b]
                triangular_adjacency(ia, ja, ib, jb) && add_edge!(g, a, b)
            end

            push!(results, (g, all_phys))
        end
    end

    @info "Generated $(length(results)) triangular subset UDGs ($(nx)×$(ny) inner grid, inner vertices $(min_inner)–$(max_k))"
    @info "pinset in generated graphs: [1,2,3,4]"
    return _process_and_save_graphs(results, path)
end

function _call_shortg(temp_path::String, mapping_file::String)
    if Sys.which("shortg") === nothing
        # Make shortg optional: log and signal caller to fallback
        @warn "Optional tool `shortg` not found; skipping canonicalization/dedup."
        return false
    else
        @info "shortg found in PATH; running canonicalization"
    end
    run(pipeline(`shortg -v -u $(temp_path)`, stderr=mapping_file))
    return true
end

function _process_and_save_graphs(results::Vector{Tuple{SimpleGraph{T}, Vector{Tuple{Float64, Float64}}}}, path::String) where T
    # If shortg is unavailable, save directly without dedup/canonicalization.
    if Sys.which("shortg") === nothing
        @warn "`shortg` not found in PATH; saving graphs without deduplication."
        save_graph(results, path)
        return path
    end

    mapping_file = tempname()
    temp_path = tempname()

    original_coords = getindex.(results, 2)
    save_graph(results, temp_path; g6_only=true)

    ok = _call_shortg(temp_path, mapping_file)
    if !ok
        # Fallback path if shortg disappeared between checks or failed early
        @warn "Skipping canonicalization due to `shortg` issue; saving raw graphs."
        try
            isfile(mapping_file) && rm(mapping_file)
            isfile(temp_path) && rm(temp_path)
        catch
        end
        save_graph(results, path)
        return path
    end

    canonical_to_original, _ = _parse_shortg_mapping(mapping_file)
    
    _write_original_representatives(temp_path, canonical_to_original, original_coords, path)

    try
        isfile(mapping_file) && rm(mapping_file)
        isfile(temp_path) && rm(temp_path)
    catch
    end
    return path
end

function _parse_shortg_mapping(filepath::String)
    canonical_to_originals = Dict{Int, Vector{Int}}()
    original_to_canonical = Dict{Int, Int}()

    for line in eachline(filepath)
        if isempty(strip(line)) || startswith(line, '>') || startswith(line, 'Z')
            continue
        end

        if occursin(":", line)
            parts = split(line, ":")
            canonical = parse(Int, strip(parts[1]))
            originals = split(strip(parts[2]))
            for orig_str in originals
                original = parse(Int, orig_str)
                push!(get!(canonical_to_originals, canonical, Int[]), original)
                original_to_canonical[original] = canonical
            end
        end
    end

    return canonical_to_originals, original_to_canonical
end

function _write_original_representatives(
    original_file::String,
    canon2orig::Dict{Int, Vector{Int}},
    orig_coords::Vector{Vector{Tuple{Float64, Float64}}},
    output_file::String
)
    original_lines = readlines(original_file)
    to_be_written = Vector{Tuple{AbstractString, Vector{Tuple{Float64, Float64}}}}()
    open(output_file, "w") do io
        for canon_line in sort(collect(keys(canon2orig)))
            orig_line = canon2orig[canon_line][1]
            push!(to_be_written, (original_lines[orig_line], orig_coords[orig_line]))
        end
    end
    save_graph(to_be_written, output_file)
end
