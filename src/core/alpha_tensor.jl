"""
    calculate_alpha_tensor(graph, boundary_vertices)

Compute the alpha tensor of a graph with given boundary vertices using tropical
tensor network contraction (Definition 3.2 in the paper).

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

Compute the reduced alpha tensor α̃(R) of a graph with given boundary vertices
(Definition 3.4 in the paper).

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

# Arguments
- `t1::AbstractArray{T}`: First reduced alpha tensor (finite values or `-Inf`)
- `t2::AbstractArray{T}`: Second reduced alpha tensor (finite values or `-Inf`)

# Returns
- `Tuple{Bool, Real}`: `(is_valid, constant_offset)`
"""
function is_diff_by_constant(t1::AbstractArray{T}, t2::AbstractArray{T}) where T <: Real
    size(t1) == size(t2) || throw(DimensionMismatch("input tensors must have the same size, got $(size(t1)) and $(size(t2))"))
    x = NaN
    for (a, b) in zip(t1, t2)
        if isinf(a) && isinf(b)
            continue
        end
        if isinf(a) || isinf(b)
            return false, 0
        end
        if isnan(x)
            x = a - b
        elseif x != a - b
            return false, 0
        end
    end
    return true, x
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
