using GadgetSearch
using Test

@testset "Core Search Functions" begin
    include("core/search.jl")
end

@testset "Unweighted Search" begin
    include("core/unweighted_search.jl")
end

@testset "Alpha Tensor" begin
    include("core/alpha_tensor.jl")
end

@testset "UDG Generation" begin
    include("graphio/udg.jl")
end

@testset "Graph Loader" begin
    include("graphio/graphloader.jl")
end

@testset "Save Graph" begin
    include("graphio/savegraph.jl")
end

@testset "Visualization" begin
    include("utils/visualize.jl")
end

@testset "Rule IO" begin
    include("utils/ruleio.jl")
end

@testset "Gadget Utilities" begin
    include("utils/gadget.jl")
end