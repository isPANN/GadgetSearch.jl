using Test
using GadgetSearch
using Graphs

@testset "save_graph - Vector{SimpleGraph}" begin
    g1 = SimpleGraph(2)
    add_edge!(g1, 1, 2)
    g2 = SimpleGraph(3)
    add_edge!(g2, 1, 2); add_edge!(g2, 2, 3)

    tmp = tempname() * ".g6"
    try
        result = save_graph([g1, g2], tmp)
        @test result == tmp
        @test isfile(tmp)
        lines = readlines(tmp)
        @test length(lines) == 2
        # Each line should be a non-empty string (g6 encoding)
        @test all(!isempty, lines)
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

@testset "save_graph - Vector{Tuple{SimpleGraph, positions}}" begin
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    pos = [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]

    tmp = tempname() * ".g6"
    try
        result = save_graph([(g, pos)], tmp)
        @test result == tmp
        @test isfile(tmp)
        line = readline(tmp)
        # Should contain coordinate data
        @test occursin("(", line)
        @test occursin(")", line)
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

@testset "save_graph - g6_only flag" begin
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    pos = [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]

    tmp = tempname() * ".g6"
    try
        save_graph([(g, pos)], tmp; g6_only=true)
        line = readline(tmp)
        # g6_only: no coordinate parentheses
        @test !occursin("(", line)
        @test !isempty(strip(line))
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

@testset "save_graph - Vector{Tuple{AbstractString, positions}}" begin
    g6 = "Bw"   # complete triangle: 3 vertices, 3 edges
    pos = [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]

    tmp = tempname() * ".g6"
    try
        result = save_graph([(g6, pos)], tmp)
        @test result == tmp
        @test isfile(tmp)
        line = readline(tmp)
        @test startswith(line, "Bw")
        @test occursin("(0.0,0.0)", line)
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

@testset "save_graph - empty vector" begin
    tmp = tempname() * ".g6"
    try
        save_graph(SimpleGraph{Int}[], tmp)
        @test isfile(tmp)
        @test isempty(read(tmp, String))
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end
