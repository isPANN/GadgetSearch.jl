function _parse_vertex_count(s::AbstractString)::Tuple{Int, Int}
    isempty(s) && throw(ArgumentError("Empty Graph6 string"))

    if s[1] != '~'
        return Int(s[1]) - 63, 2
    elseif length(s) >= 4 && s[2] != '~'
        nv = 0
        for i in 2:4
            nv = (nv << 6) + (Int(s[i]) - 63)
        end
        return nv, 5
    else
        length(s) >= 8 || throw(ArgumentError("Malformed Graph6 vertex-count prefix"))
        nv = 0
        for i in 3:8
            nv = (nv << 6) + (Int(s[i]) - 63)
        end
        return nv, 9
    end
end

function _encode_vertex_count(nv::Int)::String
    nv >= 0 || throw(ArgumentError("vertex count must be non-negative"))

    if nv <= 62
        return string(Char(nv + 63))
    elseif nv <= 258_047
        payload = Char[]
        push!(payload, '~')
        for shift in (12, 6, 0)
            push!(payload, Char(((nv >> shift) & 0x3f) + 63))
        end
        return String(payload)
    elseif nv <= 68_719_476_735
        payload = Char[]
        push!(payload, '~')
        push!(payload, '~')
        for shift in (30, 24, 18, 12, 6, 0)
            push!(payload, Char(((nv >> shift) & 0x3f) + 63))
        end
        return String(payload)
    end

    throw(ArgumentError("graph6 supports at most 68719476735 vertices, got $nv"))
end

"""
    graph_to_g6(g::SimpleGraph{Int}; include_header::Bool=false) -> String

Encode a graph to graph6 text. By default this returns the canonical graph6
payload without the `>>graph6<<` prefix so it can be used directly in
`GraphDataset`.
"""
function graph_to_g6(g::SimpleGraph{Int}; include_header::Bool=false)::String
    nv_g = nv(g)
    payload = IOBuffer()
    print(payload, _encode_vertex_count(nv_g))

    bits = Bool[]
    sizehint!(bits, nv_g * (nv_g - 1) ÷ 2)
    @inbounds for col in 2:nv_g, row in 1:(col - 1)
        push!(bits, has_edge(g, row, col))
    end

    nbits = length(bits)
    i = 1
    while i <= nbits
        value = 0
        for _ in 1:6
            value <<= 1
            if i <= nbits && bits[i]
                value |= 0x01
            end
            i += 1
        end
        print(payload, Char(value + 63))
    end

    encoded = String(take!(payload))
    if include_header
        return ">>graph6<<" * encoded
    end
    return encoded
end

function _parse_g6_string(s::AbstractString, temp_bitvec::BitVector)::SimpleGraph{Int}
    startswith(s, ">>graph6<<") && (s = s[11:end])
    isempty(s) && throw(ArgumentError("Empty Graph6 string"))

    nv, pos = _parse_vertex_count(s)
    nv <= 0 && return SimpleGraph(0)

    g = SimpleGraph(nv)

    nbits = nv * (nv - 1) ÷ 2
    if length(temp_bitvec) < nbits
        resize!(temp_bitvec, nbits)
    end
    fill!(temp_bitvec, false)

    bit_idx = 1
    @inbounds for c in s[pos:end]
        val = Int(c) - 63
        for i in 5:-1:0
            bit_idx > nbits && break
            temp_bitvec[bit_idx] = ((val >> i) & 1) == 1
            bit_idx += 1
        end
    end

    bit_idx = 1
    @inbounds for col in 2:nv, row in 1:(col - 1)
        temp_bitvec[bit_idx] && add_edge!(g, row, col)
        bit_idx += 1
    end

    return g
end
