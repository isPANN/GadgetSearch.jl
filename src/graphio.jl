function _bv2int(x::BitVector)
    @assert(length(x) <= 8 * sizeof(Int))
    acc = 0
    for i in eachindex(x)
        acc = (acc << 1) + x[i]
    end
    return acc
end

function _int2bv(n::Int, k::Int)
    bitstr = lstrip(bitstring(n), '0')
    l = length(bitstr)
    padding = k - l
    bv = falses(k)
    for i in 1:l
        bv[padding + i] = (bitstr[i] == '1')
    end
    return bv
end

function _g6_R(_x::BitVector)::Vector{UInt8}
    k = length(_x)
    padding = cld(k, 6) * 6 - k
    x = vcat(_x, falses(padding))
    nbytes = div(length(x), 6)
    bytevec = Vector{UInt8}(undef, nbytes)   # uninitialized data!
    for i in 1:nbytes
        xslice = x[((i - 1) * 6 + 1):(i * 6)]

        intslice = 0
        for bit in xslice
            intslice = (intslice << 1) + bit
        end
        intslice += 63
        bytevec[i] = intslice
    end
    return UInt8.(bytevec)
end

_g6_R(n::Int, k::Int) = _g6_R(_int2bv(n, k))

function _g6_Rp(bytevec::Vector{UInt8})
    nbytes = length(bytevec)
    x = BitVector()
    for byte in bytevec
        bits = _int2bv(byte - 63, 6)
        x = vcat(x, bits)
    end
    return x
end

function _g6_N(x::Integer)::Vector{UInt8}
    if (x < 0) || (x > 68719476735)
        error("x must satisfy 0 <= x <= 68719476735")
    elseif (x <= 62)
        nvec = [x + 63]
    elseif (x <= 258047)
        nvec = vcat([0x7e], _g6_R(x, 18))
    else
        nvec = vcat([0x7e; 0x7e], _g6_R(x, 36))
    end
    return UInt8.(nvec)
end

function _g6_Np(N::Vector{UInt8})
    if N[1] < 0x7e
        return (Int(N[1] - 63), N[2:end])
    elseif N[2] < 0x7e
        return (_bv2int(_g6_Rp(N[2:4])), N[5:end])
    else
        return (_bv2int(_g6_Rp(N[3:8])), N[9:end])
    end
end

function g6string_to_graph(s::AbstractString)::SimpleGraph
    if startswith(s, ">>graph6<<")
        s = s[11:end]
    end
    V = Vector{UInt8}(s)
    (nv, rest) = _g6_Np(V)
    bitvec = _g6_Rp(rest)
    n = 0
    g = Graphs.SimpleGraph(nv)
    for i in 2:nv, j in 1:(i - 1)
        n += 1
        if bitvec[n]
            add_edge!(g, j, i)
        end
    end
    return g
end

function read_g6_file(file_path::String)::Vector{SimpleGraph}
    results = SimpleGraph[]
    open(file_path, "r") do io
        for line in eachline(io)
            graph_data = strip(line)
            if isempty(graph_data)
                continue
            end
            push!(results, g6string_to_graph(graph_data))
        end
    end
    return results
end

