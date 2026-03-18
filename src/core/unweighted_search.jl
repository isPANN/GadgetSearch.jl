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
    equivalent_representations(graph, boundary; max_added_vertices=0, preserve_boundary_roles=true) -> Vector{Tuple{SimpleGraph{Int}, Vector{Int}}}

Generate all equivalent representations of a target graph with its boundary.

Two `(graph, boundary)` pairs are equivalent if they encode the same gadget rule.
Currently this includes:

- edge-subdivision expansions of the target graph, with up to
  `max_added_vertices` inserted path vertices distributed across the target
  edges
- optional boundary permutations that yield distinct reduced alpha tensors
  when `preserve_boundary_roles=false`

The input pair is always included as the first element of the returned vector.

By default, the input boundary ordering is preserved so boundary pin roles stay
fixed during crossing-gadget search. Set `preserve_boundary_roles=false` to
recover the older numbering-agnostic behavior.

This function belongs to the graph-theory layer and is independently useful for
inspection, testing, and future equivalence strategies such as logical flips.

# Arguments
- `graph::SimpleGraph{Int}`: The target graph R
- `boundary::Vector{Int}`: Boundary vertices of R

# Keywords
- `max_added_vertices::Int=0`: Maximum number of inserted internal vertices to
  distribute across the target edges when generating expanded representations

# Returns
- `Vector{Tuple{SimpleGraph{Int}, Vector{Int}}}`: All distinct equivalent
  `(graph, boundary)` pairs. The original `(graph, boundary)` is guaranteed to
  be the first element.
"""
function equivalent_representations(
    graph::SimpleGraph{Int},
    boundary::Vector{Int},
    ;
    max_added_vertices::Int=0,
    preserve_boundary_roles::Bool=true,
)
    max_added_vertices >= 0 || throw(ArgumentError("max_added_vertices must be non-negative"))

    seeds = Tuple{SimpleGraph{Int}, Vector{Int}}[]
    for expanded_graph in _subdivision_graph_variants(graph, max_added_vertices)
        append!(
            seeds,
            _boundary_equivalent_representations(
                expanded_graph,
                boundary;
                preserve_boundary_roles=preserve_boundary_roles,
            ),
        )
    end

    return _dedup_equivalent_representations(seeds)
end

function _subdivision_graph_variants(graph::SimpleGraph{Int}, max_added_vertices::Int)
    max_added_vertices >= 0 || throw(ArgumentError("max_added_vertices must be non-negative"))

    edge_list = sort([(min(src(e), dst(e)), max(src(e), dst(e))) for e in edges(graph)])
    allocations = _subdivision_allocations(length(edge_list), max_added_vertices)
    return [_subdivide_graph_edges(graph, edge_list, counts) for counts in allocations]
end

function _subdivision_allocations(edge_count::Int, max_added_vertices::Int)
    allocations = Vector{Vector{Int}}()
    current = zeros(Int, edge_count)

    function recurse(edge_idx::Int, remaining::Int)
        if edge_idx > edge_count
            push!(allocations, copy(current))
            return
        end

        for inserted in 0:remaining
            current[edge_idx] = inserted
            recurse(edge_idx + 1, remaining - inserted)
        end
    end

    recurse(1, max_added_vertices)
    return allocations
end

function _subdivide_graph_edges(
    graph::SimpleGraph{Int},
    edge_list::AbstractVector{<:Tuple{Int, Int}},
    inserted_counts::AbstractVector{<:Integer},
)
    length(edge_list) == length(inserted_counts) || throw(ArgumentError(
        "edge_list and inserted_counts must have the same length"))
    all(==(0), inserted_counts) && return graph

    expanded = SimpleGraph(nv(graph))

    for ((u, v), count) in zip(edge_list, inserted_counts)
        count < 0 && throw(ArgumentError("inserted edge-subdivision counts must be non-negative"))

        if count == 0
            add_edge!(expanded, u, v)
            continue
        end

        previous = u
        for _ in 1:count
            add_vertex!(expanded)
            next_vertex = nv(expanded)
            add_edge!(expanded, previous, next_vertex)
            previous = next_vertex
        end
        add_edge!(expanded, previous, v)
    end

    return expanded
end

function _boundary_equivalent_representations(
    graph::SimpleGraph{Int},
    boundary::Vector{Int},
    ;
    preserve_boundary_roles::Bool,
)
    if preserve_boundary_roles
        return [(graph, copy(boundary))]
    end

    return _boundary_permutation_representations(graph, boundary)
end

function _boundary_permutation_representations(
    graph::SimpleGraph{Int},
    boundary::Vector{Int},
)
    return [(graph, collect(perm)) for perm in Combinatorics.permutations(boundary)]
end

function _dedup_equivalent_representations(
    seeds::AbstractVector{<:Tuple{SimpleGraph{Int}, Vector{Int}}},
)
    isempty(seeds) && throw(ArgumentError("seeds must be non-empty"))

    seen = Set{Vector{Float64}}()
    result = Tuple{SimpleGraph{Int}, Vector{Int}}[]

    for (graph, boundary) in seeds
        reduced = vec(Float64.(content.(calculate_reduced_alpha_tensor(graph, boundary))))
        if reduced ∉ seen
            push!(seen, reduced)
            push!(result, (graph, boundary))
        end
    end

    return result
end

function _offset_from_pattern(
    target_graph::SimpleGraph{Int},
    pattern_graph::SimpleGraph{Int},
    pattern_boundary::Vector{Int},
    pattern_reduced::AbstractVector{<:Real},
    offset_cache::IdDict{SimpleGraph{Int}, Union{Nothing, Float64}},
)
    if haskey(offset_cache, target_graph)
        return offset_cache[target_graph]
    end

    if target_graph === pattern_graph
        offset_cache[target_graph] = 0.0
        return 0.0
    end

    target_reduced = vec(Float64.(content.(calculate_reduced_alpha_tensor(target_graph, pattern_boundary))))
    equivalent, offset = is_diff_by_constant(target_reduced, pattern_reduced)
    result = equivalent ? float(offset) : nothing
    offset_cache[target_graph] = result
    return result
end

function _build_target_data(
    representations::AbstractVector{<:Tuple{SimpleGraph{Int}, Vector{Int}}};
    include_logical_flips::Bool=true,
)
    isempty(representations) && throw(ArgumentError("representations must be non-empty"))

    k = length(representations[1][2])
    pattern_graph, pattern_boundary = representations[1]
    pattern_tensor = Float64.(content.(calculate_reduced_alpha_tensor(pattern_graph, pattern_boundary)))
    pattern_reduced = vec(pattern_tensor)
    all(isinf.(pattern_reduced)) && throw(ArgumentError(
        "target graph has an entirely -Inf reduced alpha tensor"))

    flip_patterns =
        if include_logical_flips
            generate_flip_patterns(k)
        else
            [(Int[], "no-flip")]
        end
    seen = Set{Vector{Float64}}()
    structural_offset_cache = IdDict{SimpleGraph{Int}, Union{Nothing, Float64}}()
    target_data_list = NamedTuple{(:graph, :boundary, :reduced, :mask, :offset_from_pattern, :maps_to_pattern, :flip_mask, :flip_desc), Tuple{SimpleGraph{Int}, Vector{Int}, Vector{Float64}, BigInt, Union{Nothing, Float64}, Bool, Vector{Int}, String}}[]

    for (target_graph, target_boundary) in representations
        length(target_boundary) == k || throw(ArgumentError(
            "all representations must have the same boundary size, got $k and $(length(target_boundary))"))

        target_tensor = Float64.(content.(calculate_reduced_alpha_tensor(target_graph, target_boundary)))
        all(isinf.(target_tensor)) && throw(ArgumentError(
            "target graph has an entirely -Inf reduced alpha tensor"))

        for (flip_mask, flip_desc) in flip_patterns
            transformed_tensor = apply_flip_to_tensor(target_tensor, flip_mask)
            target_reduced = vec(Float64.(transformed_tensor))
            target_reduced ∈ seen && continue

            maps_to_pattern = false
            offset_from_pattern =
                if isempty(flip_mask)
                    structural_offset = _offset_from_pattern(
                        target_graph,
                        pattern_graph,
                        pattern_boundary,
                        pattern_reduced,
                        structural_offset_cache,
                    )
                    maps_to_pattern = structural_offset !== nothing
                    structural_offset
                else
                    nothing
                end

            push!(seen, target_reduced)
            push!(target_data_list, (
                graph=target_graph,
                boundary=target_boundary,
                reduced=target_reduced,
                mask=inf_mask(target_reduced),
                offset_from_pattern=offset_from_pattern,
                maps_to_pattern=maps_to_pattern,
                flip_mask=copy(flip_mask),
                flip_desc=flip_desc,
            ))
        end
    end

    return target_data_list
end

"""
    _make_unweighted_filter(representations; prefilter)

Internal: build a filter closure from a set of equivalent target representations.

Pre-computes the reduced alpha tensor and `-Inf` mask for each representation,
then returns a closure that checks candidate graphs against all of them.
"""
function _make_unweighted_filter(
    representations::AbstractVector{<:Tuple{SimpleGraph{Int}, Vector{Int}}};
    prefilter::Bool=true,
    include_logical_flips::Bool=true,
)
    isempty(representations) && throw(ArgumentError("representations must be non-empty"))

    k = length(representations[1][2])
    pattern_graph, pattern_boundary = representations[1]
    target_data_list = _build_target_data(
        representations;
        include_logical_flips=include_logical_flips,
    )

    apply_prefilter = prefilter && !_allows_unpinned_components(representations)

    return function(candidate::SimpleGraph{Int}, pos, pin_set)
        vertex_pool = something(pin_set, 1:nv(candidate))

        if apply_prefilter && !pins_prefilter(candidate, vertex_pool)
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
                    total_offset =
                        if td.maps_to_pattern
                            float(constant_offset + something(td.offset_from_pattern))
                        else
                            float(constant_offset)
                        end
                    return UnweightedGadget(
                        pattern_graph,
                        candidate,
                        collect(boundary),
                        total_offset,
                        pos,
                    )
                end
            end
        end

        return nothing
    end
end

function _allows_unpinned_components(
    representations::AbstractVector{<:Tuple{SimpleGraph{Int}, Vector{Int}}},
)
    return any(representations) do (graph, boundary)
        boundary_set = Set(boundary)
        any(component -> !any(in(boundary_set), component), connected_components(graph))
    end
end

"""
    search_unweighted_gadgets(target_graph, target_boundary, loader; max_added_vertices, preserve_boundary_roles, prefilter, limit, max_results)

Search for unweighted gadget replacements of `target_graph` by iterating over a
`GraphLoader`.

Internally, [`equivalent_representations`](@ref) is called to expand the target
into all distinct search representations, including edge-subdivision expansions
and logical-flip tensor variants. By default the boundary ordering is preserved,
so pin `1/2/3/4` keeps its canonical role during crossing search. The
subdivision budget is controlled explicitly by `max_added_vertices`, so the
target semantics do not depend on the searched loader.

# Arguments
- `target_graph::SimpleGraph{Int}`: The pattern graph R
- `target_boundary::Vector{Int}`: Boundary vertices of R
- `loader::GraphLoader`: Graph dataset to search over

# Keywords
- `max_added_vertices::Int=0`: Maximum number of inserted internal vertices to
  distribute across target edges when generating equivalent representations
- `preserve_boundary_roles::Bool=true`: Keep the input boundary ordering fixed
  when generating equivalent representations. Set to `false` to allow distinct
  boundary permutations during matching.
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
    max_added_vertices::Int=0,
    preserve_boundary_roles::Bool=true,
    prefilter::Bool=true,
    include_logical_flips::Bool=true,
    limit::Union{Int,Nothing}=nothing,
    max_results::Union{Int,Nothing}=nothing
)
    max_added_vertices >= 0 || throw(ArgumentError("max_added_vertices must be non-negative"))

    total = limit === nothing ? length(loader) : min(length(loader), limit)
    reprs = equivalent_representations(
        target_graph,
        target_boundary;
        max_added_vertices=max_added_vertices,
        preserve_boundary_roles=preserve_boundary_roles,
    )
    filter_fn = _make_unweighted_filter(
        reprs;
        prefilter=prefilter,
        include_logical_flips=include_logical_flips,
    )
    results = UnweightedGadget[]

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


