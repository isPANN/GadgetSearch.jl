using Test
using GadgetSearch
using Graphs
using JSON3

@testset "save_graph - Vector{SimpleGraph} JSONL" begin
    g1 = SimpleGraph(2)
    add_edge!(g1, 1, 2)
    g2 = SimpleGraph(3)
    add_edge!(g2, 1, 2); add_edge!(g2, 2, 3)

    tmp = tempname() * ".jsonl"
    try
        result = save_graph([g1, g2], tmp)
        @test result == tmp
        @test isfile(tmp)
        lines = readlines(tmp)
        @test length(lines) == 2
        obj1 = JSON3.read(lines[1])
        @test haskey(obj1, :g6)
        @test !haskey(obj1, :pos)
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

@testset "save_graph - Vector{Tuple{SimpleGraph, positions}} JSONL" begin
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    pos = [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]

    tmp = tempname() * ".jsonl"
    try
        result = save_graph([(g, pos)], tmp)
        @test result == tmp
        lines = readlines(tmp)
        obj = JSON3.read(lines[1])
        @test haskey(obj, :g6)
        @test haskey(obj, :pos)
        @test length(obj[:pos]) == 3
        @test obj[:pos][1] == [0.0, 0.0]
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

@testset "save_graph - Vector{Tuple{AbstractString, positions}} JSONL" begin
    g6 = "Bw"
    pos = [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]

    tmp = tempname() * ".jsonl"
    try
        result = save_graph([(g6, pos)], tmp)
        @test result == tmp
        lines = readlines(tmp)
        obj = JSON3.read(lines[1])
        @test obj[:g6] == "Bw"
        @test obj[:pos][1] == [0.0, 0.0]
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

@testset "save_graph - empty vector" begin
    tmp = tempname() * ".jsonl"
    try
        save_graph(SimpleGraph{Int}[], tmp)
        @test isfile(tmp)
        @test isempty(read(tmp, String))
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end
