"""
    UnweightedGadget

Result of an unweighted gadget search (Definition 3.6 / Theorem 3.7 in the paper).

# Fields
- `pattern_graph::SimpleGraph{Int}`: The target graph R
- `replacement_graph::SimpleGraph{Int}`: The candidate graph R' that replaces R
- `boundary_vertices::Vector{Int}`: Boundary vertices of `replacement_graph`
- `constant_offset::Float64`: Constant offset between the reduced alpha tensors of `pattern_graph` and `replacement_graph`
- `target_index::Int`: Index of the matched target graph, defaults to `1` for single-target search
- `pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}`: Vertex positions of `replacement_graph`
"""
struct UnweightedGadget
    pattern_graph::SimpleGraph{Int}
    replacement_graph::SimpleGraph{Int}
    boundary_vertices::Vector{Int}
    constant_offset::Float64
    target_index::Int
    pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}
end

function UnweightedGadget(
    pattern_graph::SimpleGraph{Int},
    replacement_graph::SimpleGraph{Int},
    boundary_vertices::Vector{Int},
    constant_offset::Float64,
    pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}
)
    return UnweightedGadget(pattern_graph, replacement_graph, boundary_vertices, constant_offset, 1, pos)
end

"""
    make_unweighted_filter(target_graph, target_boundary)

Create a filter closure for unweighted gadget search.

For a single target, this is a thin wrapper around the multi-target method.

# Arguments
- `target_graph::SimpleGraph{Int}`: The pattern graph R
- `target_boundary::Vector{Int}`: Boundary vertices of R

# Keywords
- `prefilter::Bool=false`: Whether to reject candidates whose allowed pin set
  misses an entire connected component before tensor computation

# Returns
- `Function`: Closure `(candidate, pos, pin_set) -> UnweightedGadget | nothing`
"""
function make_unweighted_filter(
    target_graph::SimpleGraph{Int},
    target_boundary::Vector{Int};
    prefilter::Bool=false
)
    return make_unweighted_filter([(target_graph, target_boundary)]; prefilter=prefilter)
end

"""
    make_unweighted_filter(targets; prefilter=true)

Create a filter closure for unweighted gadget search over multiple targets.

The reduced alpha tensors and `-Inf` masks of all targets are pre-computed once.
Candidate boundaries are then compared only against targets with matching
boundary size, using `inf_mask` for fast rejection before the constant-difference
check.

# Arguments
- `targets::AbstractVector{<:Tuple{SimpleGraph{Int}, Vector{Int}}}`: Target
  graphs paired with their boundary vertices

# Keywords
- `prefilter::Bool=true`: Whether to reject candidates whose allowed pin set
  misses an entire connected component before tensor computation

# Returns
- `Function`: Closure `(candidate, pos, pin_set) -> UnweightedGadget | nothing`
"""
function make_unweighted_filter(
    targets::AbstractVector{<:Tuple{SimpleGraph{Int}, Vector{Int}}};
    prefilter::Bool=true
)
    isempty(targets) && throw(ArgumentError("targets must be non-empty"))

    target_groups = Dict{Int, Vector{NamedTuple{(:pattern_graph, :target_index, :reduced, :mask), Tuple{SimpleGraph{Int}, Int, Vector{Float64}, BigInt}}}}()
    group_order = Int[]

    for (target_index, (target_graph, target_boundary)) in enumerate(targets)
        target_reduced = vec(Float64.(content.(calculate_reduced_alpha_tensor(target_graph, target_boundary))))
        all(isinf.(target_reduced)) && throw(ArgumentError("target graph at index $target_index has an entirely -Inf reduced alpha tensor"))

        k = length(target_boundary)
        if !haskey(target_groups, k)
            target_groups[k] = NamedTuple{(:pattern_graph, :target_index, :reduced, :mask), Tuple{SimpleGraph{Int}, Int, Vector{Float64}, BigInt}}[]
            push!(group_order, k)
        end

        push!(target_groups[k], (
            pattern_graph=target_graph,
            target_index=target_index,
            reduced=target_reduced,
            mask=inf_mask(target_reduced),
        ))
    end

    return function(candidate::SimpleGraph{Int}, pos, pin_set)
        vertex_pool = something(pin_set, 1:nv(candidate))

        if prefilter && !pins_prefilter(candidate, vertex_pool)
            return nothing
        end

        for k in group_order
            (nv(candidate) < k || length(vertex_pool) < k) && continue

            for boundary in Combinatorics.combinations(vertex_pool, k)
                candidate_reduced = vec(Float64.(content.(calculate_reduced_alpha_tensor(candidate, boundary))))
                all(isinf.(candidate_reduced)) && continue

                candidate_mask = inf_mask(candidate_reduced)
                for target_data in target_groups[k]
                    candidate_mask == target_data.mask || continue

                    valid, constant_offset = is_diff_by_constant(candidate_reduced, target_data.reduced)
                    if valid
                        return UnweightedGadget(
                            target_data.pattern_graph,
                            candidate,
                            collect(boundary),
                            float(constant_offset),
                            target_data.target_index,
                            pos,
                        )
                    end
                end
            end
        end

        return nothing
    end
end

"""
    search_unweighted_gadgets(target_graph, target_boundary, loader; limit, max_results)

Search for unweighted gadget replacements of a single target graph.

This is a thin wrapper around the multi-target method and returns the same
`Vector{UnweightedGadget}` result type.

# Arguments
- `target_graph::SimpleGraph{Int}`: The pattern graph R
- `target_boundary::Vector{Int}`: Boundary vertices of R
- `loader::GraphLoader`: Graph dataset to search over

# Keywords
- `prefilter::Bool=false`: Whether to reject candidates whose allowed pin set
  misses an entire connected component before tensor computation
- `limit::Union{Int,Nothing}=nothing`: Maximum number of graphs to examine
- `max_results::Union{Int,Nothing}=nothing`: Stop after finding this many results

# Returns
- `Vector{UnweightedGadget}`
"""
function search_unweighted_gadgets(
    target_graph::SimpleGraph{Int},
    target_boundary::Vector{Int},
    loader::GraphLoader;
    prefilter::Bool=false,
    limit::Union{Int,Nothing}=nothing,
    max_results::Union{Int,Nothing}=nothing
)
    return search_unweighted_gadgets(
        [(target_graph, target_boundary)],
        loader;
        prefilter=prefilter,
        limit=limit,
        max_results=max_results,
    )
end

"""
    search_unweighted_gadgets(targets, loader; prefilter=true, limit=nothing, max_results=nothing)

Search for unweighted gadget replacements of multiple target graphs by iterating
over a `GraphLoader`.

For each candidate graph, boundary combinations are tested against all targets
with the same boundary size. Matches are returned as `UnweightedGadget`s, and
`target_index` records which target matched.

# Arguments
- `targets::AbstractVector{<:Tuple{SimpleGraph{Int}, Vector{Int}}}`: Target
  graphs paired with their boundary vertices
- `loader::GraphLoader`: Graph dataset to search over

# Keywords
- `prefilter::Bool=true`: Whether to reject candidates whose allowed pin set
  misses an entire connected component before tensor computation
- `limit::Union{Int,Nothing}=nothing`: Maximum number of graphs to examine
- `max_results::Union{Int,Nothing}=nothing`: Stop after finding this many results

# Returns
- `Vector{UnweightedGadget}`
"""
function search_unweighted_gadgets(
    targets::AbstractVector{<:Tuple{SimpleGraph{Int}, Vector{Int}}},
    loader::GraphLoader;
    prefilter::Bool=true,
    limit::Union{Int,Nothing}=nothing,
    max_results::Union{Int,Nothing}=nothing
)
    filter_fn = make_unweighted_filter(targets; prefilter=prefilter)
    results = UnweightedGadget[]
    total = limit === nothing ? length(loader) : min(length(loader), limit)

    @showprogress for key in Iterators.take(keys(loader), total)
        result = filter_fn(loader[key], loader.layout[key], loader.pinset)
        result === nothing && continue
        push!(results, result)
        max_results !== nothing && length(results) >= max_results && break
    end
    return results
end

"""
    calculate_alpha_tensor(graph, boundary_vertices)

Compute the alpha tensor of a graph with given boundary vertices using tropical
tensor network contraction.

The alpha tensor `alpha(R)` is a rank-`|boundary|` tensor where each element
`alpha(R)_s` is the size of the largest independent set of `R` with boundary
vertices fixed to configuration `s`, or `-Inf` if `s` violates the independent
set constraint.

# Arguments
- `graph::SimpleGraph{Int}`: The graph R
- `boundary_vertices::Vector{Int}`: The boundary vertex indices (1-indexed)

# Returns
- `Array{<:Tropical}`: Tropical tensor of shape `(2, 2, ..., 2)` with `length(boundary_vertices)` dimensions
"""
function calculate_alpha_tensor(graph::SimpleGraph{Int}, boundary_vertices::Vector{Int})
    return solve(GenericTensorNetwork(IndependentSet(graph), openvertices=boundary_vertices), SizeMax())
end

"""
    calculate_reduced_alpha_tensor(graph, boundary_vertices)

Compute the reduced alpha tensor α̃(R) of a graph with given boundary vertices.

The reduced alpha tensor is obtained by applying `mis_compactify!` to the alpha tensor:
a boundary configuration `s` is set to `-Inf` if there exists a subset configuration
`s' ⊆ s` (fewer boundary vertices selected) that achieves an equal or better MIS size.
This removes dominated configurations, leaving only the "essential" boundary behaviors.

# Arguments
- `graph::SimpleGraph{Int}`: The graph R
- `boundary_vertices::Vector{Int}`: The boundary vertex indices (1-indexed)

# Returns
- `Array{<:Tropical}`: Reduced tropical tensor of shape `(2, 2, ..., 2)`
"""
function calculate_reduced_alpha_tensor(graph::SimpleGraph{Int}, boundary_vertices::Vector{Int})
    return mis_compactify!(calculate_alpha_tensor(graph, boundary_vertices))
end

"""
    inf_mask(tensor)

Return a bitmask encoding the positions of `-Inf` entries in `tensor`.

The first linearized tensor entry corresponds to the least-significant bit.
"""
function inf_mask(tensor::AbstractArray)
    mask = BigInt(0)
    for (i, value) in enumerate(tensor)
        entry = value isa Tropical ? content(value) : value
        if isinf(entry) && entry < 0
            mask |= BigInt(1) << (i - 1)
        end
    end
    return mask
end

"""
    pins_prefilter(g, pins)

Return `true` when every connected component of `g` contains at least one pin.
"""
function pins_prefilter(g::SimpleGraph{Int}, pins::AbstractVector{<:Integer})
    isempty(pins) && return false

    n = nv(g)
    unique_pins = unique(Int.(pins))
    length(unique_pins) == length(pins) || throw(ArgumentError("pins must be unique"))
    all(1 .<= unique_pins .<= n) || throw(ArgumentError("pins must be valid vertex indices for a graph with $n vertices"))
    pinset = Set(unique_pins)
    return all(component -> any(in(pinset), component), connected_components(g))
end

"""
    is_diff_by_constant(t1, t2)

Check whether two reduced alpha tensors differ by a constant (Theorem 3.7 in the paper).

Two gadgets R and R' are valid MIS replacements if and only if their reduced alpha tensors
differ only by a constant offset. This function tests that condition.

Returns `(true, c)` if `t1[i] - t2[i] == c` for all finite entries, and both tensors have
`-Inf` at exactly the same positions. Returns `(false, 0)` otherwise.

Throws `ArgumentError` if both tensors are entirely `-Inf`.

# Arguments
- `t1::AbstractArray{T}`: First reduced alpha tensor (finite values or `-Inf`)
- `t2::AbstractArray{T}`: Second reduced alpha tensor (finite values or `-Inf`)

# Returns
- `Tuple{Bool, Real}`: `(is_valid, constant_offset)`
"""
function is_diff_by_constant(t1::AbstractArray{T}, t2::AbstractArray{T}) where T <: Real
    size(t1) == size(t2) || throw(DimensionMismatch("input tensors must have the same size, got $(size(t1)) and $(size(t2))"))
    any(isinf.(t1) .⊻ isinf.(t2)) && return false, zero(T)
    d = filter(isfinite, t1 .- t2)
    isempty(d) && throw(ArgumentError("input tensors must contain at least one finite entry"))
    return all(==(first(d)), d), first(d)
end

"""
    is_gadget_replacement(g1, g2, open_vertices1, open_vertices2)

Check whether gadget `g2` (with boundary `open_vertices2`) is a valid MIS replacement
for gadget `g1` (with boundary `open_vertices1`), i.e., their reduced alpha tensors
differ only by a constant (Theorem 3.7 in the paper).

# Returns
- `Tuple{Bool, Real}`: `(is_valid, constant_offset)` where `constant_offset = α̃(g2) - α̃(g1)`
"""
function is_gadget_replacement(g1::SimpleGraph{Int}, g2::SimpleGraph{Int},
                                open_vertices1::Vector{Int}, open_vertices2::Vector{Int})
    t1 = content.(calculate_reduced_alpha_tensor(g1, open_vertices1))
    t2 = content.(calculate_reduced_alpha_tensor(g2, open_vertices2))
    return is_diff_by_constant(t2, t1)
end


