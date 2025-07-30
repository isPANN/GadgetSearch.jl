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

