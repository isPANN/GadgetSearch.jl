# ============================================================================
# Unweighted Gadget Types
# ============================================================================

"""
    UnweightedGadget

Result of an unweighted gadget search.
Stores the pattern graph R, replacement graph R', boundary vertices,
constant offset between reduced alpha tensors, and optional vertex positions.
"""
struct UnweightedGadget
    pattern_graph::SimpleGraph{Int}
    replacement_graph::SimpleGraph{Int}
    boundary_vertices::Vector{Int}
    constant_offset::Float64
    pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}
end

# ============================================================================
# Target Data
# ============================================================================

function _build_target_data(
    pattern_graph::SimpleGraph{Int},
    pattern_boundary::Vector{Int};
    include_logical_flips::Bool=true,
)
    k = length(pattern_boundary)
    pattern_tensor = Float64.(content.(calculate_reduced_alpha_tensor(pattern_graph, pattern_boundary)))
    pattern_reduced = vec(pattern_tensor)
    all(isinf.(pattern_reduced)) && error(
        "target graph has an entirely -Inf reduced alpha tensor")
    flip_patterns =
        if include_logical_flips
            generate_flip_patterns(k)
        else
            [(Int[], "no-flip")]
        end
    seen = Set{Vector{Float64}}()
    target_data_list = NamedTuple{(:reduced, :mask, :flip_mask, :flip_desc), Tuple{Vector{Float64}, BigInt, Vector{Int}, String}}[]
    for (flip_mask, flip_desc) in flip_patterns
        transformed_tensor = apply_flip_to_tensor(pattern_tensor, flip_mask)
        target_reduced = vec(Float64.(transformed_tensor))
        target_reduced ∈ seen && continue
        push!(seen, target_reduced)
        push!(target_data_list, (
            reduced=target_reduced,
            mask=inf_mask(target_reduced),
            flip_mask=copy(flip_mask),
            flip_desc=flip_desc,
        ))
    end
    return target_data_list
end

# ============================================================================
# Unweighted Filter Construction
# ============================================================================

"""
    _make_unweighted_filter(pattern_graph, pattern_boundary; prefilter, include_logical_flips)

Build a filter closure that checks candidate graphs against the target pattern.
"""
function _make_unweighted_filter(
    pattern_graph::SimpleGraph{Int},
    pattern_boundary::Vector{Int};
    prefilter::Bool=true,
    include_logical_flips::Bool=true,
)
    k = length(pattern_boundary)
    target_data_list = _build_target_data(
        pattern_graph, pattern_boundary;
        include_logical_flips=include_logical_flips,
    )
    apply_prefilter = prefilter && all(
        component -> any(in(Set(pattern_boundary)), component),
        Graphs.connected_components(pattern_graph),
    )
    return function(candidate::SimpleGraph{Int}, pos, pin_set)
        vertex_pool = something(pin_set, 1:Graphs.nv(candidate))
        if apply_prefilter && !pins_prefilter(candidate, vertex_pool)
            return nothing
        end
        (Graphs.nv(candidate) < k || length(vertex_pool) < k) && return nothing
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

# ============================================================================
# Unweighted Search
# ============================================================================

"""
    search_unweighted_gadgets(target_graph, target_boundary, loader; kwargs...)

Search for unweighted gadget replacements of `target_graph` by iterating over a `GraphLoader`.
"""
function search_unweighted_gadgets(
    target_graph::SimpleGraph{Int},
    target_boundary::Vector{Int},
    loader::GraphLoader;
    prefilter::Bool=true,
    include_logical_flips::Bool=true,
    limit::Union{Int,Nothing}=nothing,
    max_results::Union{Int,Nothing}=nothing,
)
    total = limit === nothing ? length(loader) : min(length(loader), limit)
    filter_fn = _make_unweighted_filter(
        target_graph, target_boundary;
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

# ============================================================================
# Alpha Tensor Functions
# ============================================================================

"""
    calculate_alpha_tensor(graph, boundary_vertices) -> Array{<:Tropical}

Compute the alpha tensor α(R) via tropical tensor network contraction.
Element `α(R)_s` is the maximum independent set size with boundary fixed to `s`,
or `-Inf` if `s` violates the independent set constraint.
"""
function calculate_alpha_tensor(graph::SimpleGraph{Int}, boundary_vertices::Vector{Int})
    return solve(GenericTensorNetwork(IndependentSet(graph), openvertices=boundary_vertices), SizeMax())
end

"""
    calculate_reduced_alpha_tensor(graph, boundary_vertices) -> Array{<:Tropical}

Compute the reduced alpha tensor α̃(R) by applying `mis_compactify!` to `α(R)`.
Dominated boundary configurations (where a subset achieves equal or better MIS size)
are set to `-Inf`.
"""
function calculate_reduced_alpha_tensor(graph::SimpleGraph{Int}, boundary_vertices::Vector{Int})
    return mis_compactify!(calculate_alpha_tensor(graph, boundary_vertices))
end

# ============================================================================
# Tensor Utility Functions
# ============================================================================

"""
    inf_mask(tensor) -> BigInt

Bitmask encoding `-Inf` positions in `tensor` (LSB = first linear index).
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
    n = Graphs.nv(g)
    unique_pins = unique(Int.(pins))
    length(unique_pins) == length(pins) || error("pins must be unique")
    all(1 .<= unique_pins .<= n) || error("pins must be valid vertex indices for a graph with $n vertices")
    pinset = Set(unique_pins)
    return all(component -> any(in(pinset), component), Graphs.connected_components(g))
end

"""
    is_diff_by_constant(t1, t2) -> (Bool, Real)

Return `(true, c)` if `t1 - t2 == c` at all finite entries and both tensors share
the same `-Inf` pattern; `(false, 0)` otherwise.
"""
function is_diff_by_constant(t1::AbstractArray{T}, t2::AbstractArray{T}) where T <: Real
    size(t1) == size(t2) || error("input tensors must have the same size, got $(size(t1)) and $(size(t2))")
    any(isinf.(t1) .⊻ isinf.(t2)) && return false, zero(T)
    d = filter(isfinite, t1 .- t2)
    isempty(d) && error("input tensors must contain at least one finite entry")
    return all(==(first(d)), d), first(d)
end

"""
    is_gadget_replacement(g1, g2, open_vertices1, open_vertices2) -> (Bool, Real)

Check whether `g2` is a valid MIS replacement for `g1`, i.e., their reduced alpha
tensors differ only by a constant. Returns `(is_valid, α̃(g2) - α̃(g1))`.
"""
function is_gadget_replacement(g1::SimpleGraph{Int}, g2::SimpleGraph{Int},
                                open_vertices1::Vector{Int}, open_vertices2::Vector{Int})
    t1 = content.(calculate_reduced_alpha_tensor(g1, open_vertices1))
    t2 = content.(calculate_reduced_alpha_tensor(g2, open_vertices2))
    return is_diff_by_constant(t2, t1)
end
