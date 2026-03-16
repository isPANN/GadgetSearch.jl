using Test
using GadgetSearch
using Graphs
using JSON3

# Helper: build a simple Rydberg gadget
function _make_rydberg_gadget()
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    constraint = GadgetSearch.TruthTableConstraint(BitMatrix([1 0; 0 1]))
    pos = [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]
    return GadgetSearch.Gadget(GadgetSearch.RydbergModel, constraint, g, [1, 2], [2.0, 2.0, 1.0], pos)
end

@testset "check_gadget - basic Rydberg" begin
    gadget = _make_rydberg_gadget()
    # _return_info=true returns a String, not nothing
    info = check_gadget(gadget; _return_info=true)
    @test info isa String
    @test occursin("Rydberg", info)
    @test occursin("Max energy", info)
    @test occursin("Ground states", info)
end

@testset "check_gadget_rydberg and check_gadget_qubo" begin
    gadget = _make_rydberg_gadget()
    rydberg_info = check_gadget_rydberg(gadget; _return_info=true)
    qubo_info = check_gadget_qubo(gadget; _return_info=true)
    @test rydberg_info isa String
    @test qubo_info isa String
    @test occursin("Rydberg", rydberg_info)
    @test occursin("QUBO", qubo_info)
end

@testset "check_gadget - default (logs, returns nothing)" begin
    gadget = _make_rydberg_gadget()
    result = check_gadget(gadget)
    @test result === nothing
end

@testset "check_gadget - argument validation" begin
    g = SimpleGraph(2)
    add_edge!(g, 1, 2)
    constraint = GadgetSearch.TruthTableConstraint(BitMatrix([1 0; 0 1]))
    # Wrong number of weights
    bad_gadget = GadgetSearch.Gadget(GadgetSearch.RydbergModel, constraint, g, [1, 2], [1.0], nothing)
    @test_throws ArgumentError check_gadget(bad_gadget; _return_info=true)

    # Pin out of range
    bad_pins = GadgetSearch.Gadget(GadgetSearch.RydbergModel, constraint, g, [1, 5], [1.0, 1.0], nothing)
    @test_throws ArgumentError check_gadget(bad_pins; _return_info=true)
end

@testset "save_results_to_json - Rydberg gadget" begin
    gadget = _make_rydberg_gadget()
    tmp = tempname() * ".json"
    try
        result = save_results_to_json([gadget], tmp)
        @test result == tmp
        @test isfile(tmp)
        data = JSON3.read(read(tmp, String))
        @test length(data) == 1
        entry = data[1]
        @test haskey(entry, "ground_states")
        @test haskey(entry, "pins")
        @test haskey(entry, "graph")
        @test entry["pins"] == [1, 2]
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

@testset "save_results_to_json - includes positions" begin
    gadget = _make_rydberg_gadget()
    tmp = tempname() * ".json"
    try
        save_results_to_json([gadget], tmp)
        data = JSON3.read(read(tmp, String))
        @test haskey(data[1]["graph"], "positions")
        @test length(data[1]["graph"]["positions"]) == 3
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

@testset "save_results_to_json - no positions" begin
    g = SimpleGraph(2)
    add_edge!(g, 1, 2)
    constraint = GadgetSearch.TruthTableConstraint(BitMatrix([1 0; 0 1]))
    gadget_no_pos = GadgetSearch.Gadget(GadgetSearch.RydbergModel, constraint, g, [1, 2], [1.0, 1.0], nothing)

    tmp = tempname() * ".json"
    try
        save_results_to_json([gadget_no_pos], tmp)
        data = JSON3.read(read(tmp, String))
        @test !haskey(data[1]["graph"], "positions")
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

@testset "save_results_to_json - empty list" begin
    tmp = tempname() * ".json"
    try
        save_results_to_json(GadgetSearch.Gadget[], tmp)
        data = JSON3.read(read(tmp, String))
        @test length(data) == 0
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

@testset "save_results_to_json - QUBO gadget includes edge weights" begin
    g = SimpleGraph(3)
    add_edge!(g, 1, 2); add_edge!(g, 2, 3)
    constraint = GadgetSearch.TruthTableConstraint(BitMatrix([1 0; 0 1]))
    edge_list = [(1, 2), (2, 3)]
    gadget = GadgetSearch.Gadget(GadgetSearch.QUBOModel, constraint, g, [1, 2],
                                  [1.0, 1.0, 1.0], [0.5, 0.5], edge_list, nothing)

    tmp = tempname() * ".json"
    try
        save_results_to_json([gadget], tmp)
        data = JSON3.read(read(tmp, String))
        @test haskey(data[1], "edge_weights")
        @test length(data[1]["edge_weights"]) == 2
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end
