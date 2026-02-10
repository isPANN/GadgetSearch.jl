using GadgetSearch
using Test

@testset "Search Strategy" begin
    include("core/search.jl")
end

@testset "UDG Generation" begin
    include("graphio/udg.jl")
end

@testset "Core Search Functions" begin
    include("core/search.jl")
end

@testset "RydbergUnweightedModel" begin
    include("core/search_unweighted.jl")
end

@testset "Visualization" begin
    include("utils/visualize.jl")
end