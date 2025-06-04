using GadgetSearch
using Test

# @testset "utils.jl" begin
#     include("utils.jl")
# end

# @testset "generateUDG.jl" begin
#     include("generateUDG.jl")
# end

# @testset "graphsearch.jl" begin
#     include("graphsearch.jl")
# end

# @testset "dataloaders.jl" begin
#     include("dataloaders.jl")
# end
# @testset "Gadget Search" begin
#     include("graphio.jl")
# end

@testset "Search Strategy" begin
    include("search.jl")
end

@testset "GraphLoader" begin
    include("graphloader.jl")
end