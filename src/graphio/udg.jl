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
get_radius(::Triangular) = 1.5

get_physical_positions(::Square, pos::Vector{Tuple{Int, Int}}) = Vector{Tuple{Float64, Float64}}(pos)
function get_physical_positions(::Triangular, pos::Vector{Tuple{Int, Int}})
    return Vector{Tuple{Float64, Float64}}([(Float64(x) - Float64(y)/2, Float64(y) * sqrt(3)/2) for (x, y) in pos])
end

# Grid shape types — describe graph connectivity on a lattice
abstract type GridShape end
struct GridSG <: GridShape end    # Grid Sub Graph — 4 neighbors (no diagonals)
struct KSG <: GridShape end       # King Sub Graph — 8 neighbors (all diagonals)
struct TLSG <: GridShape end      # Triangular Lattice Sub Graph — 6 neighbors

function parse_shape(s::String)
    s == "grid" && return GridSG()
    s == "KSG" && return KSG()
    s == "TLSG" && return TLSG()
    error("Unknown shape: $s")
end

shape_name(::GridSG) = "grid"
shape_name(::KSG) = "KSG"
shape_name(::TLSG) = "TLSG"

get_shape(::Square) = KSG()
get_shape(::Triangular) = TLSG()

function get_pin_positions(::LatticeType, nx::Int, ny::Int)
    top    = Tuple{Int, Int}[(1, y)     for y in 2:ny+1]
    bottom = Tuple{Int, Int}[(nx+2, y)  for y in 2:ny+1]
    left   = Tuple{Int, Int}[(x, 1)     for x in 2:nx+1]
    right  = Tuple{Int, Int}[(x, ny+2)  for x in 2:nx+1]
    return top, bottom, left, right
end

function get_inner_points(::LatticeType, nx::Int, ny::Int)
    return vec(Tuple{Int, Int}[(x, y) for x in 2:nx+1, y in 2:ny+1])
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
    generate_full_grid_graph(lattice::LatticeType, nx::Int, ny::Int; path::String="grid.jsonl") -> String

Generate complete graphs on a grid lattice (without boundary expansion) and save to file.
All vertices are connected to each other, regardless of distance.

# Arguments
- `lattice`: Type of lattice (Square or Triangular) - determines physical positions
- `nx`: Number of grid points in x direction
- `ny`: Number of grid points in y direction
- `path`: Output file path for saving graphs (default: "grid.jsonl")

# Returns
- `String`: Path to the saved graph file

# Details
Generates a single complete graph on the nx×ny grid. The physical positions
follow the lattice geometry, but edges connect all pairs of vertices.
"""
function generate_full_grid_graph(lattice::LatticeType, nx::Int, ny::Int; path::String="grid.jsonl")
    int_pos = vec(Tuple{Int, Int}[(x, y) for x in 1:nx, y in 1:ny])
    n = length(int_pos)
    g = complete_graph(n)
    sname = shape_name(get_shape(lattice))
    results = Tuple{SimpleGraph{Int}, String, Vector{Tuple{Int, Int}}}[(g, sname, int_pos)]
    save_graph(results, path)
    @info "Generated complete graph with $n vertices on $(nx)×$(ny) grid"
    return path
end

"""
    generate_full_grid_udg(lattice::LatticeType, nx::Int, ny::Int; path::String="udg.jsonl") -> String

Generate unit disk graphs on a grid lattice with four boundary pins and save to file.

# Arguments
- `lattice`: Type of lattice (Square or Triangular)
- `nx`: Number of inner positions in x direction
- `ny`: Number of inner positions in y direction
- `path`: Output file path for saving graphs (default: "udg.jsonl")

# Returns
- `String`: Path to the saved graph file

# Details
Generates all possible UDGs by placing pins on boundary positions and 
connecting vertices within unit distance on the specified lattice type.
"""
function generate_full_grid_udg(lattice::LatticeType, nx::Int, ny::Int; path::String="udg.jsonl")
    top_candidates, bottom_candidates, left_candidates, right_candidates = get_pin_positions(lattice, nx, ny)
    inner_points = get_inner_points(lattice, nx, ny)
    radius = get_radius(lattice)
    sname = shape_name(get_shape(lattice))
    results = Tuple{SimpleGraph{Int}, String, Vector{Tuple{Int, Int}}}[]
    for top in top_candidates, bottom in bottom_candidates,
        left in left_candidates, right in right_candidates
        int_pos = vcat([top, right, bottom, left], inner_points)
        physical = get_physical_positions(lattice, int_pos)
        g = unit_disk_graph(physical, radius)
        push!(results, (g, sname, int_pos))
    end
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

function _process_and_save_graphs(results::Vector{Tuple{SimpleGraph{T}, String, Vector{Tuple{Int, Int}}}}, path::String) where T
    if Sys.which("shortg") === nothing
        @warn "`shortg` not found in PATH; saving graphs without deduplication."
        save_graph(results, path)
        return path
    end
    mapping_file = tempname()
    temp_jsonl = tempname() * ".jsonl"
    temp_g6 = tempname() * ".g6"
    original_data = [(r[2], r[3]) for r in results]
    save_graph(results, temp_jsonl)
    export_g6(temp_jsonl, temp_g6)
    ok = _call_shortg(temp_g6, mapping_file)
    if !ok
        @warn "Skipping canonicalization due to `shortg` issue; saving raw graphs."
        try
            isfile(mapping_file) && rm(mapping_file)
            isfile(temp_jsonl) && rm(temp_jsonl)
            isfile(temp_g6) && rm(temp_g6)
        catch
        end
        save_graph(results, path)
        return path
    end
    canonical_to_original, _ = _parse_shortg_mapping(mapping_file)
    _write_original_representatives(temp_jsonl, canonical_to_original, original_data, path)
    try
        isfile(mapping_file) && rm(mapping_file)
        isfile(temp_jsonl) && rm(temp_jsonl)
        isfile(temp_g6) && rm(temp_g6)
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
    orig_data::Vector{Tuple{String, Vector{Tuple{Int, Int}}}},
    output_file::String
)
    original_lines = readlines(original_file)
    to_be_written = Vector{Tuple{AbstractString, String, Vector{Tuple{Int, Int}}}}()
    for canon_line in sort(collect(keys(canon2orig)))
        orig_line = canon2orig[canon_line][1]
        obj = JSON3.read(original_lines[orig_line])
        shape, int_pos = orig_data[orig_line]
        push!(to_be_written, (String(obj[:g6]), shape, int_pos))
    end
    save_graph(to_be_written, output_file)
end
