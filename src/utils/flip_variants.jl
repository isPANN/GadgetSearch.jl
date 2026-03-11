"""
    make_flip_aware_multi_target_filter(base_targets, flip_patterns; prefilter=true, permute_pins=true)

Create a filter that checks candidates against base targets AND their flip variants.

# Arguments
- `base_targets`: Vector of (graph, boundary, description) tuples
- `flip_patterns`: Vector of (flip_mask, flip_desc) tuples

# Keywords
- `prefilter`: Whether to apply the fast pin connectivity prefilter
- `permute_pins`: Whether to also compare against all permutations of target pin order

# Returns
- Filter function and target descriptions vector
"""
function make_flip_aware_multi_target_filter(base_targets, flip_patterns;
                                             prefilter=true,
                                             permute_pins=true)
    # Pre-compute all target tensors (base × flip combinations)
    target_data = @NamedTuple{graph::SimpleGraph{Int}, boundary::Vector{Int}, reduced::Array{Float64}, mask::UInt}[]
    target_descs = String[]

    for (g, b, desc) in base_targets
        base_tensor = content.(calculate_reduced_alpha_tensor(g, b))
        all(isinf.(base_tensor)) && continue

        perms = permute_pins ? collect(Combinatorics.permutations(1:length(b))) : [collect(1:length(b))]

        for (flip_mask, flip_desc) in flip_patterns
            # Apply flip to tensor
            flipped_tensor = apply_flip_to_tensor(base_tensor, flip_mask)

            for perm in perms
                reduced = perm == collect(1:length(b)) ? flipped_tensor : permutedims(flipped_tensor, Tuple(perm))
                push!(target_data, (graph=g, boundary=b, reduced=reduced, mask=inf_mask(reduced)))

                if permute_pins
                    push!(target_descs, "$desc-$flip_desc-perm$(join(perm))")
                else
                    push!(target_descs, "$desc-$flip_desc")
                end
            end
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
            valid && return (target_index=i, gadget=UnweightedGadget(td.graph, candidate, pins, float(offset), pos))
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

