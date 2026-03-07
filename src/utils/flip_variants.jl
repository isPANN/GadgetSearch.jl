"""
Utilities for generating flip variants of target graphs for crossing gadget search.
"""

using GenericTensorNetworks: content

"""
    generate_flip_patterns()

Generate non-equivalent flip patterns for 4-pin CROSS graph (considering symmetry).

Returns vector of (flip_mask, description) tuples.
"""
function generate_flip_patterns()
    return [
        (Int[], "no-flip"),
        ([1], "flip-pin1"),
        ([1,2], "flip-pin1-2"),
        ([1,3], "flip-pin1-3"),
        ([1,2,3,4], "flip-all")
    ]
end

"""
    generate_extended_cross()

Generate CROSS variants with inserted nodes on edges.

Returns vector of (graph, boundary, description) tuples.
"""
function generate_extended_cross()
    variants = Tuple{SimpleGraph{Int}, Vector{Int}, String}[]

    # Base CROSS
    cross = SimpleGraph(4)
    add_edge!(cross, 1, 3)
    add_edge!(cross, 2, 4)
    push!(variants, (cross, [1,2,3,4], "base"))

    # Variant 1: Insert node 5 on edge 1-3
    g1 = SimpleGraph(5)
    add_edge!(g1, 1, 5); add_edge!(g1, 5, 3)
    add_edge!(g1, 2, 4)
    push!(variants, (g1, [1,2,3,4], "ext5v-13"))

    # Variant 2: Insert node 5 on edge 2-4
    g2 = SimpleGraph(5)
    add_edge!(g2, 1, 3)
    add_edge!(g2, 2, 5); add_edge!(g2, 5, 4)
    push!(variants, (g2, [1,2,3,4], "ext5v-24"))

    # Variant 3: Insert nodes on both edges
    g3 = SimpleGraph(6)
    add_edge!(g3, 1, 5); add_edge!(g3, 5, 3)
    add_edge!(g3, 2, 6); add_edge!(g3, 6, 4)
    push!(variants, (g3, [1,2,3,4], "ext6v-both"))

    return variants
end

"""
    make_flip_aware_multi_target_filter(base_targets, flip_patterns; prefilter=true)

Create a filter that checks candidates against base targets AND their flip variants.

# Arguments
- `base_targets`: Vector of (graph, boundary, description) tuples
- `flip_patterns`: Vector of (flip_mask, flip_desc) tuples

# Returns
- Filter function and target descriptions vector
"""
function make_flip_aware_multi_target_filter(base_targets, flip_patterns; prefilter=true)
    # Pre-compute all target tensors (base × flip combinations)
    target_data = []
    target_descs = String[]

    for (g, b, desc) in base_targets
        base_tensor = content.(calculate_reduced_alpha_tensor(g, b))
        all(isinf.(base_tensor)) && continue

        for (flip_mask, flip_desc) in flip_patterns
            # Apply flip to tensor
            flipped_tensor = apply_flip_to_tensor(base_tensor, flip_mask)
            push!(target_data, (graph=g, boundary=b, reduced=flipped_tensor, mask=inf_mask(flipped_tensor)))
            push!(target_descs, "$desc-$flip_desc")
        end
    end

    k = length(base_targets[1][2])

    filter_fn = function(candidate::SimpleGraph{Int}, pos, pin_set)
        nv(candidate) < k && return nothing
        pins = something(pin_set, collect(1:nv(candidate)))

        if prefilter
            pins_prefilter(candidate, pins) || return nothing
        end

        candidate_reduced = content.(calculate_reduced_alpha_tensor(candidate, pins))
        all(isinf.(candidate_reduced)) && return nothing
        candidate_mask = inf_mask(candidate_reduced)

        for (i, td) in enumerate(target_data)
            candidate_mask != td.mask && continue
            valid, offset = is_diff_by_constant(candidate_reduced, td.reduced)
            valid && return MultiTargetResult(i, UnweightedGadget(td.graph, candidate, pins, float(offset), pos))
        end
        return nothing
    end

    return filter_fn, target_descs
end

"""
    apply_flip_to_tensor(tensor::Array{T,N}, pins_to_flip::Vector{Int}) where {T,N}

Apply logical flip to specified dimensions of a tensor.
"""
function apply_flip_to_tensor(tensor::Array{T,N}, pins_to_flip::Vector{Int}) where {T,N}
    isempty(pins_to_flip) && return tensor

    result = similar(tensor)
    dims = size(tensor)

    for idx in CartesianIndices(tensor)
        new_idx = Tuple(idx)
        for pin in pins_to_flip
            new_idx = Base.setindex(new_idx, dims[pin] + 1 - new_idx[pin], pin)
        end
        result[idx] = tensor[new_idx...]
    end

    return result
end

