"""
Utilities for generating logical-flip target variants for crossing gadget search.
"""

"""
    generate_flip_patterns(pin_count::Int=4)

Generate boundary-pin flip patterns for a target with `pin_count` boundary pins.

For the 4-pin crossing target we reuse the symmetry-reduced set explored in the
historical crossing-gadget branch. For other pin counts we currently keep only
the no-flip case to avoid exploding the target search space unintentionally.
"""
function generate_flip_patterns(pin_count::Int=4)
    pin_count < 0 && throw(ArgumentError("pin_count must be non-negative"))

    if pin_count == 4
        return [
            (Int[], "no-flip"),
            ([1], "flip-pin1"),
            ([1, 2], "flip-pin1-2"),
            ([1, 3], "flip-pin1-3"),
            ([1, 2, 3, 4], "flip-all"),
        ]
    end

    return [(Int[], "no-flip")]
end

"""
    apply_flip_to_tensor(tensor, pins_to_flip)

Apply a logical flip to the specified tensor dimensions. Flipping pin `j`
reverses the `j`-th dimension, which corresponds to swapping boundary values 0
and 1 for that pin.
"""
function apply_flip_to_tensor(
    tensor::AbstractArray{T, N},
    pins_to_flip::AbstractVector{<:Integer},
) where {T, N}
    isempty(pins_to_flip) && return copy(tensor)

    result = similar(tensor)
    dims = size(tensor)

    for idx in CartesianIndices(tensor)
        source_idx = Tuple(idx)
        flipped_idx = source_idx
        for pin in pins_to_flip
            1 <= pin <= N || throw(ArgumentError("flip pin index $pin is out of bounds for a rank-$N tensor"))
            flipped_idx = Base.setindex(flipped_idx, dims[pin] + 1 - flipped_idx[pin], pin)
        end
        result[idx] = tensor[flipped_idx...]
    end

    return result
end
