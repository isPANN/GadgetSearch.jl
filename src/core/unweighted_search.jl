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
    make_unweighted_filter(target_graph, target_boundary)

Create a filter closure for unweighted gadget search.

Pre-computes `α̃(target_graph)` once, then returns a closure that checks each
candidate graph by trying all boundary vertex combinations of the same size.

# Arguments
- `target_graph::SimpleGraph{Int}`: The pattern graph R
- `target_boundary::Vector{Int}`: Boundary vertices of R

# Returns
- `Function`: Closure `(candidate, pos, pin_set) -> UnweightedGadget | nothing`
"""
function make_unweighted_filter(target_graph::SimpleGraph{Int}, target_boundary::Vector{Int})
    target_reduced = content.(calculate_reduced_alpha_tensor(target_graph, target_boundary))
    all(isinf.(target_reduced)) && throw(ArgumentError("target graph has an entirely -Inf reduced alpha tensor"))
    k = length(target_boundary)

    return function(candidate::SimpleGraph{Int}, pos, pin_set)
        nv(candidate) < k && return nothing
        vertex_pool = something(pin_set, 1:nv(candidate))
        for boundary in Combinatorics.combinations(vertex_pool, k)
            candidate_reduced = content.(calculate_reduced_alpha_tensor(candidate, boundary))
            all(isinf.(candidate_reduced)) && continue
            valid, constant_offset = is_diff_by_constant(candidate_reduced, target_reduced)
            valid && return UnweightedGadget(target_graph, candidate, boundary, float(constant_offset), pos)
        end
        return nothing
    end
end

"""
    search_unweighted_gadgets(target_graph, target_boundary, loader; limit, max_results)

Search for unweighted gadget replacements of `target_graph` by iterating over a `GraphLoader`.

For each candidate graph, tries all boundary vertex combinations of size
`length(target_boundary)` and checks if the reduced alpha tensors differ by a constant
(Theorem 3.7). Returns on the first valid boundary found per candidate.

# Arguments
- `target_graph::SimpleGraph{Int}`: The pattern graph R
- `target_boundary::Vector{Int}`: Boundary vertices of R
- `loader::GraphLoader`: Graph dataset to search over

# Keywords
- `limit::Union{Int,Nothing}=nothing`: Maximum number of graphs to examine
- `max_results::Union{Int,Nothing}=nothing`: Stop after finding this many results

# Returns
- `Vector{UnweightedGadget}`
"""
function search_unweighted_gadgets(
    target_graph::SimpleGraph{Int},
    target_boundary::Vector{Int},
    loader::GraphLoader;
    limit::Union{Int,Nothing}=nothing,
    max_results::Union{Int,Nothing}=nothing
)
    filter_fn = make_unweighted_filter(target_graph, target_boundary)
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
    is_diff_by_constant(t1, t2)

Check whether two reduced alpha tensors differ by a constant (Theorem 3.7 in the paper).

Two gadgets R and R' are valid MIS replacements if and only if their reduced alpha tensors
differ only by a constant offset. This function tests that condition.

Returns `(true, c)` if `t1[i] - t2[i] == c` for all finite entries, and both tensors have
`-Inf` at exactly the same positions. Returns `(false, 0)` otherwise.

Callers must ensure that at least one tensor has a finite entry (i.e., not all-Inf).

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
    @assert !isempty(d)
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

# ──────────────────────────────────────────────────────────────────────────────
# Multi-target search with pre-filters
# ──────────────────────────────────────────────────────────────────────────────

"""
    inf_mask(tensor::AbstractArray{<:Real}) -> UInt

Compute a bitmask indicating which entries of `tensor` are `-Inf`.
Bit `i-1` (0-indexed) is set if `tensor[i]` is `-Inf`.

This provides an O(1) fingerprint comparison: two reduced alpha tensors
can only differ by a constant if they have identical `-Inf` patterns.
"""
function inf_mask(tensor::AbstractArray{T}) where T <: Real
    mask = UInt(0)
    for (i, v) in enumerate(tensor)
        if isinf(v)
            mask |= UInt(1) << (i - 1)
        end
    end
    return mask
end

"""
    pins_prefilter(g::SimpleGraph, pins::Vector{Int}) -> Bool

Fast structural pre-filter (O(V+E)) that rejects candidate graphs
which cannot possibly match any crossing-type target.

Checks:
1. All boundary pins have degree ≥ 1
2. All boundary pins are in the same connected component
"""
function pins_prefilter(g::SimpleGraph, pins::Vector{Int})
    # All pins must have at least one neighbor
    for p in pins
        degree(g, p) == 0 && return false
    end
    # All pins must be in the same connected component
    for i in 2:length(pins)
        has_path(g, pins[1], pins[i]) || return false
    end
    return true
end

"""
    MultiTargetResult

Result of a multi-target unweighted gadget search.

# Fields
- `target_index::Int`: Which target graph matched (1-indexed into the targets list)
- `gadget::UnweightedGadget`: The matching gadget
"""
struct MultiTargetResult
    target_index::Int
    gadget::UnweightedGadget
end

"""
    make_multi_target_filter(targets; prefilter=true)

Create a filter closure that checks a candidate graph against **multiple**
target graphs in a single pass. The candidate's reduced alpha tensor is
computed only once, then compared against all targets via:

1. **Connectivity pre-filter** (O(V+E)): reject if pins are disconnected
2. **Inf-mask fingerprint** (O(1)): reject targets with non-matching -Inf patterns
3. **Full constant-offset check** (O(2^k)): only for fingerprint-matched targets

# Arguments
- `targets::Vector{Tuple{SimpleGraph{Int}, Vector{Int}}}`:
  List of `(pattern_graph, boundary_vertices)` pairs

# Keywords
- `prefilter::Bool=true`: Whether to apply the connectivity pre-filter

# Returns
- `Function`: Closure `(candidate, pos, pin_set) -> MultiTargetResult | nothing`
"""
function make_multi_target_filter(targets::Vector{Tuple{SimpleGraph{Int}, Vector{Int}}};
                                   prefilter::Bool=true)
    # Pre-compute reduced alpha tensors and inf-masks for all targets
    target_data = map(targets) do (g, b)
        reduced = content.(calculate_reduced_alpha_tensor(g, b))
        all(isinf.(reduced)) && throw(ArgumentError("target graph has an entirely -Inf reduced alpha tensor"))
        (graph=g, boundary=b, reduced=reduced, mask=inf_mask(reduced))
    end
    k = length(first(targets)[2])

    return function(candidate::SimpleGraph{Int}, pos, pin_set)
        nv(candidate) < k && return nothing
        pins = something(pin_set, collect(1:nv(candidate)))

        # ── Stage 1: Connectivity pre-filter (O(V+E)) ──────────────────
        if prefilter
            pins_prefilter(candidate, pins) || return nothing
        end

        # ── Stage 2: Compute candidate α̃ ONCE ──────────────────────────
        candidate_reduced = content.(calculate_reduced_alpha_tensor(candidate, pins))
        all(isinf.(candidate_reduced)) && return nothing
        candidate_mask = inf_mask(candidate_reduced)

        # ── Stage 3: Compare against all targets ────────────────────────
        for (i, td) in enumerate(target_data)
            # Fast inf-mask fingerprint check (O(1))
            candidate_mask != td.mask && continue
            # Full constant-offset check
            valid, constant_offset = is_diff_by_constant(candidate_reduced, td.reduced)
            valid && return MultiTargetResult(i,
                UnweightedGadget(td.graph, candidate, pins, float(constant_offset), pos))
        end
        return nothing
    end
end

"""
    search_multi_target_gadgets(targets, loader; prefilter, limit, max_results)

Search for unweighted gadget replacements against **multiple** target graphs
simultaneously. Each candidate graph's α̃ tensor is computed once and compared
against all targets.

# Arguments
- `targets::Vector{Tuple{SimpleGraph{Int}, Vector{Int}}}`:
  List of `(pattern_graph, boundary_vertices)` pairs
- `loader::GraphLoader`: Graph dataset to search over

# Keywords
- `prefilter::Bool=true`: Apply connectivity pre-filter before tensor computation
- `limit::Union{Int,Nothing}=nothing`: Maximum number of graphs to examine
- `max_results::Union{Int,Nothing}=nothing`: Stop after finding this many results

# Returns
- `Vector{MultiTargetResult}`
"""
function search_multi_target_gadgets(
    targets::Vector{Tuple{SimpleGraph{Int}, Vector{Int}}},
    loader::GraphLoader;
    prefilter::Bool=true,
    limit::Union{Int,Nothing}=nothing,
    max_results::Union{Int,Nothing}=nothing
)
    filter_fn = make_multi_target_filter(targets; prefilter=prefilter)
    results = MultiTargetResult[]
    total = limit === nothing ? length(loader) : min(length(loader), limit)

    @showprogress for key in Iterators.take(keys(loader), total)
        result = filter_fn(loader[key], loader.layout[key], loader.pinset)
        result === nothing && continue
        push!(results, result)
        max_results !== nothing && length(results) >= max_results && break
    end
    return results
end


