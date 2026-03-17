"""
    UnweightedGadget

Result of an unweighted gadget search (Definition 3.6 / Theorem 3.7 in the paper).

# Fields
- `pattern_graph::SimpleGraph{Int}`: The target graph R
- `replacement_graph::SimpleGraph{Int}`: The candidate graph R' that replaces R
- `boundary_vertices::Vector{Int}`: Boundary vertices of `replacement_graph`
- `constant_offset::Float64`: Constant offset between the reduced alpha tensors of `pattern_graph` and `replacement_graph`
- `pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}`: Vertex positions of `replacement_graph`
"""
struct UnweightedGadget
    pattern_graph::SimpleGraph{Int}
    replacement_graph::SimpleGraph{Int}
    boundary_vertices::Vector{Int}
    constant_offset::Float64
    pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}
end

"""
    equivalent_representations(graph, boundary) -> Vector{Tuple{SimpleGraph{Int}, Vector{Int}}}

Generate all equivalent representations of a target graph with its boundary.

Two `(graph, boundary)` pairs are equivalent if they encode the same gadget rule
— i.e. the same graph with a different boundary vertex ordering that yields a
distinct reduced alpha tensor.  The input pair is always included as the first
element of the returned vector.

This function belongs to the graph-theory layer and is independently useful for
inspection, testing, and future equivalence strategies (e.g. vertex additions).

# Arguments
- `graph::SimpleGraph{Int}`: The target graph R
- `boundary::Vector{Int}`: Boundary vertices of R

# Returns
- `Vector{Tuple{SimpleGraph{Int}, Vector{Int}}}`: All distinct equivalent
  `(graph, permuted_boundary)` pairs.  The original `(graph, boundary)` is
  guaranteed to be the first element.
"""
function equivalent_representations(
    graph::SimpleGraph{Int},
    boundary::Vector{Int},
)
    base_reduced = vec(Float64.(content.(calculate_reduced_alpha_tensor(graph, boundary))))
    seen = Set{Vector{Float64}}()
    push!(seen, base_reduced)

    result = Tuple{SimpleGraph{Int}, Vector{Int}}[(graph, boundary)]

    for perm in Combinatorics.permutations(boundary)
        perm == boundary && continue
        perm_reduced = vec(Float64.(content.(calculate_reduced_alpha_tensor(graph, perm))))
        if perm_reduced ∉ seen
            push!(seen, perm_reduced)
            push!(result, (graph, perm))
        end
    end

    return result
end

"""
    _make_unweighted_filter(representations; prefilter)

Internal: build a filter closure from a set of equivalent target representations.

Pre-computes the reduced alpha tensor and `-Inf` mask for each representation,
then returns a closure that checks candidate graphs against all of them.
"""
function _make_unweighted_filter(
    representations::AbstractVector{<:Tuple{SimpleGraph{Int}, Vector{Int}}};
    prefilter::Bool=true
)
    isempty(representations) && throw(ArgumentError("representations must be non-empty"))

    target_data_list = NamedTuple{(:pattern_graph, :reduced, :mask), Tuple{SimpleGraph{Int}, Vector{Float64}, BigInt}}[]
    k = length(representations[1][2])

    for (target_graph, target_boundary) in representations
        length(target_boundary) == k || throw(ArgumentError(
            "all representations must have the same boundary size, got $k and $(length(target_boundary))"))

        target_reduced = vec(Float64.(content.(calculate_reduced_alpha_tensor(target_graph, target_boundary))))
        all(isinf.(target_reduced)) && throw(ArgumentError(
            "target graph has an entirely -Inf reduced alpha tensor"))

        push!(target_data_list, (
            pattern_graph=target_graph,
            reduced=target_reduced,
            mask=inf_mask(target_reduced),
        ))
    end

    pattern_graph = representations[1][1]

    return function(candidate::SimpleGraph{Int}, pos, pin_set)
        vertex_pool = something(pin_set, 1:nv(candidate))

        if prefilter && !pins_prefilter(candidate, vertex_pool)
            return nothing
        end

        (nv(candidate) < k || length(vertex_pool) < k) && return nothing

        for boundary in Combinatorics.combinations(vertex_pool, k)
            candidate_reduced = vec(Float64.(content.(calculate_reduced_alpha_tensor(candidate, boundary))))
            all(isinf.(candidate_reduced)) && continue

            candidate_mask = inf_mask(candidate_reduced)
            for td in target_data_list
                candidate_mask == td.mask || continue

                valid, constant_offset = is_diff_by_constant(candidate_reduced, td.reduced)
                if valid
                    return UnweightedGadget(
                        pattern_graph,
                        candidate,
                        collect(boundary),
                        float(constant_offset),
                        pos,
                    )
                end
            end
        end

        return nothing
    end
end

"""
    search_unweighted_gadgets(target_graph, target_boundary, loader; prefilter, limit, max_results)

Search for unweighted gadget replacements of `target_graph` by iterating over a
`GraphLoader`.

Internally, [`equivalent_representations`](@ref) is called to expand the target
into all distinct boundary orderings, and candidates are matched against the
full set.  The caller sees a single-target API.

# Arguments
- `target_graph::SimpleGraph{Int}`: The pattern graph R
- `target_boundary::Vector{Int}`: Boundary vertices of R
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
    target_graph::SimpleGraph{Int},
    target_boundary::Vector{Int},
    loader::GraphLoader;
    prefilter::Bool=true,
    limit::Union{Int,Nothing}=nothing,
    max_results::Union{Int,Nothing}=nothing
)
    reprs = equivalent_representations(target_graph, target_boundary)
    filter_fn = _make_unweighted_filter(reprs; prefilter=prefilter)
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


