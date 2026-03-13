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

"""
    triangular_adjacency(i1::Int, j1::Int, i2::Int, j2::Int) -> Bool

Return whether two sites of the triangular lattice are nearest neighbors using
only integer lattice coordinates.

# Arguments
- `i1::Int`: Column index of the first site.
- `j1::Int`: Row index of the first site.
- `i2::Int`: Column index of the second site.
- `j2::Int`: Row index of the second site.

# Returns
- `Bool`: `true` if the two lattice sites share an edge in the triangular lattice.
"""
function triangular_adjacency(i1::Int, j1::Int, i2::Int, j2::Int)
    j1 == j2 && return abs(i1 - i2) == 1
    abs(j1 - j2) == 1 || return false

    if j1 > j2
        return triangular_adjacency(i2, j2, i1, j1)
    end

    return isodd(j1) ? (i2 == i1 || i2 == i1 + 1) : (i2 == i1 || i2 == i1 - 1)
end

"""
    triangular_lattice_graph(nx::Int, ny::Int) -> SimpleGraph{Int}

Build the full nearest-neighbor triangular lattice graph on an `nx × ny` grid
of lattice sites.

# Arguments
- `nx::Int`: Number of lattice columns.
- `ny::Int`: Number of lattice rows.

# Returns
- `SimpleGraph{Int}`: The graph whose vertices are lattice sites and whose
  edges connect nearest neighbors.
"""
function triangular_lattice_graph(nx::Int, ny::Int)
    if nx <= 0 || ny <= 0
        return SimpleGraph(0)
    end

    vertex_index(i, j) = i + (j - 1) * nx
    g = SimpleGraph(nx * ny)

    for j = 1:ny, i = 1:nx
        v = vertex_index(i, j)

        if i < nx
            add_edge!(g, v, vertex_index(i + 1, j))
        end

        if j < ny
            add_edge!(g, v, vertex_index(i, j + 1))

            if isodd(j) && i < nx
                add_edge!(g, v, vertex_index(i + 1, j + 1))
            elseif iseven(j) && i > 1
                add_edge!(g, v, vertex_index(i - 1, j + 1))
            end
        end
    end

    return g
end

"""
    get_radius(lattice::LatticeType) -> Float64

Return the unit-disk radius used for the chosen lattice geometry.

# Arguments
- `lattice::LatticeType`: Lattice family whose physical spacing determines the
  unit-disk threshold.

# Returns
- `Float64`: The interaction radius used when constructing unit-disk graphs.
"""
get_radius(::Square) = 1.5
get_radius(::Triangular) = 1.1

_triangular_grid_coordinates(nx::Int, ny::Int) = vec(Tuple{Int, Int}[(i, j) for i in 1:nx, j in 1:ny])

"""
    dedup_inner_subsets(inner_grid::SimpleGraph{Int}, k::Integer; use_shortg::Bool=true) -> Vector{Vector{Int}}

Enumerate all `k`-vertex subsets of `inner_grid` and, when possible, retain
one representative from each isomorphism class using `shortg`.

# Arguments
- `inner_grid::SimpleGraph{Int}`: The lattice graph whose vertex subsets are enumerated.
- `k::Integer`: Size of each subset to generate.
- `use_shortg::Bool=true`: Whether to use `shortg` for isomorphism-based deduplication
  when the executable is available in `PATH`.

# Returns
- `Vector{Vector{Int}}`: A collection of vertex subsets, each represented by the
  original vertex indices from `inner_grid`.
"""
function dedup_inner_subsets(inner_grid::SimpleGraph{Int}, k::Integer; use_shortg::Bool=true)
    0 <= k <= nv(inner_grid) || throw(ArgumentError("subset size k must satisfy 0 <= k <= nv(inner_grid)"))

    subsets = [collect(subset) for subset in Combinatorics.combinations(vertices(inner_grid), k)]
    if isempty(subsets) || !use_shortg || Sys.which("shortg") === nothing
        return subsets
    end

    subgraphs = SimpleGraph{Int}[]
    for subset in subsets
        subgraph, _ = Graphs.induced_subgraph(inner_grid, subset)
        push!(subgraphs, subgraph)
    end

    mapping_file = tempname()
    temp_path = tempname()

    try
        save_graph(subgraphs, temp_path)
        ok = _call_shortg(temp_path, mapping_file)
        ok || return subsets

        canonical_to_original, _ = _parse_shortg_mapping(mapping_file)
        return [subsets[canonical_to_original[idx][1]] for idx in sort!(collect(keys(canonical_to_original)))]
    finally
        isfile(mapping_file) && rm(mapping_file)
        isfile(temp_path) && rm(temp_path)
    end
end

get_physical_positions(::Square, pos::Vector{Tuple{Int, Int}}) = Vector{Tuple{Float64, Float64}}(pos)
function get_physical_positions(::Triangular, pos::Vector{Tuple{Int, Int}})
    h = sqrt(3) / 2
    return Vector{Tuple{Float64, Float64}}([(x + (isodd(y) ? 0.5 : 0.0), y * h) for (x, y) in pos])
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
    generate_triangular_udg_subsets(nx::Int, ny::Int;
                                    subset_sizes=0:(nx * ny),
                                    deduplicate::Bool=true,
                                    use_shortg::Bool=true,
                                    path::String="triangular_udg_subsets.g6") -> String

Generate all triangular-lattice unit-disk graphs obtained by selecting subsets
of the `nx × ny` inner lattice sites and save them to disk.

# Arguments
- `nx::Int`: Number of lattice columns in the inner triangular grid.
- `ny::Int`: Number of lattice rows in the inner triangular grid.
- `subset_sizes=0:(nx * ny)`: Iterable of subset sizes to enumerate.
- `deduplicate::Bool=true`: Whether to collapse isomorphic inner subsets before saving.
- `use_shortg::Bool=true`: Whether `shortg` may be used for deduplication when available.
- `path::String="triangular_udg_subsets.g6"`: Output file path.

# Returns
- `String`: The path where the generated dataset was written.
"""
function generate_triangular_udg_subsets(
    nx::Int,
    ny::Int;
    subset_sizes=0:(nx * ny),
    deduplicate::Bool=true,
    use_shortg::Bool=true,
    path::String="triangular_udg_subsets.g6",
)
    inner_grid = triangular_lattice_graph(nx, ny)
    grid_coords = _triangular_grid_coordinates(nx, ny)
    physical_positions = get_physical_positions(Triangular(), grid_coords)
    results = Tuple{SimpleGraph{Int}, Vector{Tuple{Float64, Float64}}}[]

    for k in subset_sizes
        0 <= k <= nv(inner_grid) || throw(ArgumentError("subset sizes must satisfy 0 <= k <= nx * ny"))

        subsets = deduplicate ?
            dedup_inner_subsets(inner_grid, k; use_shortg=use_shortg) :
            [collect(subset) for subset in Combinatorics.combinations(vertices(inner_grid), k)]

        for subset in subsets
            subset_positions = physical_positions[subset]
            subset_graph = unit_disk_graph(subset_positions, get_radius(Triangular()))
            push!(results, (subset_graph, subset_positions))
        end
    end

    save_graph(results, path)
    return path
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
