using Revise
using Graphs
using GadgetSearch

# g = SimpleGraph(14)
# add_edge!(g, 1, 2)
# add_edge!(g, 1,5)
# add_edge!(g, 3, 1)
# add_edge!(g, 2,4)
# add_edge!(g, 2,7)
# add_edge!(g, 2,5)
# add_edge!(g, 3,5)
# add_edge!(g, 3,8)
# add_edge!(g, 3,6)
# add_edge!(g, 4, 7)
# add_edge!(g, 5,7)
# add_edge!(g, 8,5)
# add_edge!(g, 8,6)
# add_edge!(g, 9,7)
# add_edge!(g, 9,10)

# add_edge!(g, 1,12)
# add_edge!(g, 4,11)
# add_edge!(g, 6,13)
# add_edge!(g, 10,14)
# # add_edge!(g, 6, 11)
# add_edge!(g, 6, 12)
# add_edge!(g, 7,12)
# add_edge!(g, 8,12)
# add_edge!(g, 9,10)
# add_edge!(g, 10,11)
# add_edge!(g, 11,12)

# add_edge!(g, 9,13)
# add_edge!(g, 9,14)
# add_edge!(g, 10,14)
# add_edge!(g, 11,14)
# add_edge!(g, 13,14)
# add_edge!(g, 13,15)
# add_edge!(g, 15,16)
# add_edge!(g, 16,17)
# add_edge!(g, 17,18)

# og= SimpleGraph(13)
# add_edge!(og, 1,2)
# add_edge!(og, 2,3)
# add_edge!(og, 3,4)
# add_edge!(og, 4,5)
# add_edge!(og, 5,6)

# add_edge!(og, 7,8)
# add_edge!(og, 8,9)
# add_edge!(og, 9,10)
# add_edge!(og, 10,11)
# add_edge!(og, 11,12)
# add_edge!(og, 12,13)
function simplegraph(edgelist::AbstractVector{Tuple{Int,Int}})
    nv = maximum(x->max(x...), edgelist)
    g = SimpleGraph(nv)
    for (i,j) in edgelist
        add_edge!(g, i, j)
    end
    return g
end

function triangular_physical_position(loc, parity=false)
    i, j = loc
    y = j * (√3 / 2)
    if parity
        x = i + (iseven(j) ? 0.5 : 0.0)
    else
        x = i + (isodd(j) ? 0.5 : 0.0)
    end
    return (x, y)
end

function triangular_unitdisk_graph(locs::AbstractVector, unit::Real, parity::Bool=false)
    # parity==false: crossing nodes on odd columns.
    n = length(locs)
    g = SimpleGraph(n)
    physical_locs = triangular_physical_position.(locs, parity)
    for i=1:n, j=i+1:n
        if sum(abs2, physical_locs[i] .- physical_locs[j]) < unit ^ 2
            add_edge!(g, i, j)
        end
    end
    return g
end

g = simplegraph([(2,1),(3,2),(5,4),(6,5)])

g = triangular_unitdisk_graph([(1,2), (2,1), (2,2), (2,3),(3,3)], 1.1, true)

weights = [2,1,3,2,1]

# og_weights = [1,2,2,2,2,1,1,2,2,2,2,2,1]
# struct Gadget{T<:Real}
#     ground_states::BitMatrix
#     graph::SimpleGraph{Int}
#     pins::Vector{Int}
#     weights::Vector{T}
#     pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}
# end
gadget = Gadget(BitMatrix([1 0 0 1;]), g, [2,1,5,4], weights, nothing)
# og_gadget = Gadget(BitMatrix([0 1 0 0;]), og, [1,7,6,13], og_weights, nothing)

GadgetSearch.check_gadget(gadget)
# GadgetSearch.check_gadget(og_gadget)

GadgetSearch.plot_graph(g,"test.pdf"; pos=nothing)