using GadgetSearch
using Test
using GLPK

@testset "utils.jl" begin
    include("utils.jl")
end

@testset "generateUDG.jl" begin
    include("generateUDG.jl")
end

@testset "graphsearch.jl" begin
    include("graphsearch.jl")
end

@testset "dataloaders.jl" begin
    include("dataloaders.jl")
end