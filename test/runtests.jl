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

@testset "Visualization" begin
    include("utils/visualize.jl")
end

@testset "Alpha Tensor" begin
    include("core/alpha_tensor.jl")
end