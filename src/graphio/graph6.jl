function _parse_vertex_count(s::AbstractString)::Tuple{Int, Int}
    if s[1] <= '~'
        return Int(s[1]) - 63, 2
    elseif length(s) > 1 && s[2] <= '~'
        nv = 0
        for i in 2:4
            nv = (nv << 6) + (Int(s[i]) - 63)
        end
        return nv, 5
    else
        nv = 0
        for i in 3:8
            nv = (nv << 6) + (Int(s[i]) - 63)
        end
        return nv, 9
    end
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