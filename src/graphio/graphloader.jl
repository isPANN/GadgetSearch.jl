function _compute_physical_positions(shape::String, grid_pos::Vector{Tuple{Int, Int}})
    if shape == "TLSG"
        return Tuple{Float64, Float64}[(Float64(x) - Float64(y)/2, Float64(y) * sqrt(3)/2) for (x, y) in grid_pos]
    else  # "KSG" or "grid"
        return Tuple{Float64, Float64}[(Float64(x), Float64(y)) for (x, y) in grid_pos]
    end
end

struct GraphDataset
    g6codes::Vector{<:AbstractString}
    shapes::Vector{Union{Nothing, String}}
    grid_positions::Vector{Union{Nothing, Vector{Tuple{Int, Int}}}}
    layouts::Vector{Union{Nothing, Vector{Tuple{Float64, Float64}}}}
    n::Int  # number of graphs
end

mutable struct GraphLoader
    # Graph dataset and parsing
    dataset::GraphDataset
    pinset::Union{Nothing, Vector{Int}}
    temp_bitvec::BitVector

    # Caching parsed graphs (optional)
    enable_cache::Bool
    max_cached::Int
    cachepath::Union{Nothing, String}
    parsed::Dict{String, SimpleGraph{Int}}
    parse_order::Vector{String}
end

struct LayoutAccessor
    cds::GraphLoader
end

function GraphDataset(g6codes::Vector{<:AbstractString},
                      shapes::Vector{Union{Nothing, String}},
                      grid_positions::Vector{Union{Nothing, Vector{Tuple{Int, Int}}}})
    n = length(g6codes)
    length(shapes) == n || throw(ArgumentError("g6codes and shapes must have the same length"))
    length(grid_positions) == n || throw(ArgumentError("g6codes and grid_positions must have the same length"))
    codes_string = String[string(code) for code in g6codes]
    layouts = Union{Nothing, Vector{Tuple{Float64, Float64}}}[]
    for i in 1:n
        if shapes[i] !== nothing && grid_positions[i] !== nothing
            push!(layouts, _compute_physical_positions(shapes[i], grid_positions[i]))
        else
            push!(layouts, nothing)
        end
    end
    return GraphDataset(codes_string, shapes, grid_positions, layouts, n)
end

function GraphDataset(g6codes::Vector{<:AbstractString})
    n = length(g6codes)
    codes_string = String[string(code) for code in g6codes]
    shapes = Union{Nothing, String}[nothing for _ in 1:n]
    grid_positions = Union{Nothing, Vector{Tuple{Int, Int}}}[nothing for _ in 1:n]
    layouts = Union{Nothing, Vector{Tuple{Float64, Float64}}}[nothing for _ in 1:n]
    return GraphDataset(codes_string, shapes, grid_positions, layouts, n)
end

function GraphDataset(path::String)
    codes = String[]
    shapes = Union{Nothing, String}[]
    grid_positions = Union{Nothing, Vector{Tuple{Int, Int}}}[]
    layouts = Union{Nothing, Vector{Tuple{Float64, Float64}}}[]
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        obj = JSON3.read(stripped)
        push!(codes, String(obj[:g6]))
        if haskey(obj, :shape) && obj[:shape] !== nothing
            shape = String(obj[:shape])
            push!(shapes, shape)
            if haskey(obj, :pos) && obj[:pos] !== nothing
                gpos = Tuple{Int, Int}[(Int(p[1]), Int(p[2])) for p in obj[:pos]]
                push!(grid_positions, gpos)
                push!(layouts, _compute_physical_positions(shape, gpos))
            else
                push!(grid_positions, nothing)
                push!(layouts, nothing)
            end
        else
            push!(shapes, nothing)
            push!(grid_positions, nothing)
            push!(layouts, nothing)
        end
    end
    return GraphDataset(codes, shapes, grid_positions, layouts, length(codes))
end


function GraphLoader(dataset::GraphDataset; cachepath::Union{Nothing, String}=nothing, max_cached::Int=10_000, enable_cache::Bool=false, pinset::Union{Nothing, Vector{Int}}=nothing)
    parsed_cache = Dict{String, SimpleGraph{Int}}()
    bitvec = BitVector(undef, 0)

    if cachepath !== nothing && isfile(cachepath)
        @info "Loading parsed graph cache from $cachepath"
        try
            loaded = deserialize(cachepath)
            if isa(loaded, Dict{String, SimpleGraph{Int}})
                parsed_cache = loaded
            else
                @warn "Cache file format unexpected, ignoring cache."
            end
        catch e
            @warn "Failed to load cache file: $e"
        end
    end

    return GraphLoader(dataset, pinset, bitvec, enable_cache, max_cached, cachepath, parsed_cache, String[])
end

function GraphLoader(path::String; cachepath::Union{Nothing, String}=nothing, max_cached::Int=10_000, enable_cache::Bool=false, pinset::Union{Nothing, Vector{Int}}=nothing)
    ds = GraphDataset(path)
    return GraphLoader(ds; cachepath=cachepath, max_cached=max_cached, enable_cache=enable_cache, pinset=pinset)
end


function Base.getindex(cds::GraphLoader, key::String)
    idx = tryparse(Int, key)
    if idx === nothing
        error("GraphLoader: String key \"$key\" is not a valid integer index.")
    end
    return cds[idx]
end

function Base.getindex(cds::GraphLoader, idx::Int)
    key = string(idx)
    if cds.enable_cache && haskey(cds.parsed, key)
        return cds.parsed[key]
    else
        g6 = cds.dataset.g6codes[idx]
        g = _parse_g6_string(g6, cds.temp_bitvec)

        if cds.enable_cache
            cds.parsed[key] = g
            push!(cds.parse_order, key)
            if length(cds.parsed) > cds.max_cached
                oldest = popfirst!(cds.parse_order)
                delete!(cds.parsed, oldest)
            end
        end
        return g
    end
end

Base.getindex(l::LayoutAccessor, key::String) = getlayout(l.cds, key)

function Base.getindex(l::LayoutAccessor, key::Int)
    return getlayout(l.cds, string(key))
end

function getlayout(cds::GraphLoader, key::String)
    idx = parse(Int, key)
    return cds.dataset.layouts[idx]
end

Base.keys(cds::GraphLoader) = 1:cds.dataset.n
Base.length(cds::GraphLoader) = cds.dataset.n

# Allow property-like access for layout: cds.layout[key]
Base.getproperty(cds::GraphLoader, name::Symbol) = name === :layout ? LayoutAccessor(cds) : getfield(cds, name)

function save_cache(cds::GraphLoader)
    if cds.cachepath === nothing
        @warn "No cache path provided; cannot save cache."
        return
    end
    @info "Saving parsed graph cache to $(cds.cachepath)"
    serialize(cds.cachepath, cds.parsed)
end

function Base.show(io::IO, loader::GraphLoader)
    nkeys = length(loader.dataset.g6codes)
    print(io, "GraphLoader with $nkeys graphs",
          loader.enable_cache ? ", cache enabled" : ", no cache",
          "")
end