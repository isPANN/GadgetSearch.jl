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

@testset "export_g6 - extracts pure g6 lines" begin
    g = SimpleGraph(3)
    add_edge!(g, 1, 2); add_edge!(g, 2, 3)
    pos = [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]

    src_file = tempname() * ".jsonl"
    dst_file = tempname() * ".g6"
    try
        save_graph([(g, pos)], src_file)
        export_g6(src_file, dst_file)
        lines = readlines(dst_file)
        @test length(lines) == 1
        @test !occursin("{", lines[1])
        @test !occursin("(", lines[1])
        bv = BitVector(undef, 0)
        parsed = GadgetSearch._parse_g6_string(lines[1], bv)
        @test nv(parsed) == 3
        @test ne(parsed) == 2
    finally
        isfile(src_file) && rm(src_file; force=true)
        isfile(dst_file) && rm(dst_file; force=true)
    end
end

@testset "export_g6 - topology-only input" begin
    g = SimpleGraph(2)
    add_edge!(g, 1, 2)

    src_file = tempname() * ".jsonl"
    dst_file = tempname() * ".g6"
    try
        save_graph([g], src_file)
        export_g6(src_file, dst_file)
        lines = readlines(dst_file)
        @test length(lines) == 1
        @test strip(lines[1]) == graph_to_g6(g)
    finally
        isfile(src_file) && rm(src_file; force=true)
        isfile(dst_file) && rm(dst_file; force=true)
    end
end

@testset "save_graph - with shape and int positions" begin
    g = SimpleGraph(3)
    add_edge!(g, 1, 2); add_edge!(g, 2, 3)
    shape = "TLSG"
    pos = [(0, 0), (1, 0), (0, 1)]

    tmp = tempname() * ".jsonl"
    try
        save_graph([(g, shape, pos)], tmp)
        lines = readlines(tmp)
        obj = JSON3.read(lines[1])
        @test obj[:shape] == "TLSG"
        @test obj[:pos] == [[0, 0], [1, 0], [0, 1]]
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

@testset "save_graph - g6 string with shape" begin
    g6 = "Bw"
    shape = "KSG"
    pos = [(0, 0), (1, 0), (0, 1)]

    tmp = tempname() * ".jsonl"
    try
        save_graph([(g6, shape, pos)], tmp)
        lines = readlines(tmp)
        obj = JSON3.read(lines[1])
        @test obj[:g6] == "Bw"
        @test obj[:shape] == "KSG"
        @test obj[:pos] == [[0, 0], [1, 0], [0, 1]]
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end
