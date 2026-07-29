using JSON3
using Graphs

include("../../app/GadgetEditor.jl")

@testset "unweighted mode computes the reduced alpha tensor" begin
    payload = JSON3.read("""
    {
      "model": "rydberg",
      "weight_mode": "unweighted",
      "nodes": [
        {"id": "3:3", "q": 3, "r": 3, "weight": 1},
        {"id": "4:3", "q": 4, "r": 3, "weight": 1},
        {"id": "5:3", "q": 5, "r": 3, "weight": 1}
      ],
      "pins": ["3:3", "5:3"]
    }
    """)

    result = GadgetEditor.compute_payload(payload)
    @test result.operation == "reduced_alpha_tensor"
    @test result.weight_mode == "unweighted"
    @test result.vertex_count == 3
    @test result.edge_count == 2
    @test result.boundary_count == 2
    @test getproperty.(result.tensor, :configuration) == ["00", "10", "01", "11"]
    @test getproperty.(result.tensor, :value) == [1.0, "-Inf", "-Inf", 2.0]
end

@testset "reduced alpha tensor can be computed without open vertices" begin
    payload = JSON3.read("""
    {
      "model": "rydberg",
      "weight_mode": "unweighted",
      "nodes": [{"id": "0:0", "q": 0, "r": 0, "weight": 1}],
      "pins": []
    }
    """)

    result = GadgetEditor.compute_payload(payload)
    @test result.operation == "reduced_alpha_tensor"
    @test result.tensor == [(configuration="", value=1.0)]
end

@testset "weighted mode uses vertex weights" begin
    payload = JSON3.read("""
    {
      "model": "rydberg",
      "weight_mode": "weighted",
      "nodes": [
        {"id": "3:3", "q": 3, "r": 3, "weight": 1},
        {"id": "4:3", "q": 4, "r": 3, "weight": 3},
        {"id": "5:3", "q": 5, "r": 3, "weight": 1}
      ],
      "pins": ["3:3", "5:3"]
    }
    """)

    result = GadgetEditor.compute_payload(payload)
    @test result.operation == "ground_states"
    @test result.weight_mode == "weighted"
    @test result.max_energy == 3.0
    @test result.observed == ["00"]
end

@testset "unweighted mode accepts more than 32 vertices" begin
    nodes = [
        (id="$(index):0", q=index, r=0, weight=1.0)
        for index in 1:33
    ]
    payload = (
        model="rydberg",
        weight_mode="unweighted",
        nodes=nodes,
        pins=String[],
    )

    result = GadgetEditor.compute_payload(payload)
    @test result.operation == "reduced_alpha_tensor"
    @test result.vertex_count == 33
    @test result.tensor == [(configuration="", value=17.0)]
end

@testset "weighted mode reports its 32-vertex enumeration limit" begin
    nodes = [
        (id="$(index):0", q=index, r=0, weight=1.0)
        for index in 1:33
    ]
    payload = (
        model="rydberg",
        weight_mode="weighted",
        nodes=nodes,
        pins=String[],
    )

    @test_throws ArgumentError GadgetEditor.gadget_from_payload(payload)
end

@testset "weighted mode requires positive weights" begin
    payload = (
        model="rydberg",
        weight_mode="weighted",
        nodes=[(id="0:0", q=0, r=0, weight=0.0)],
        pins=String[],
    )

    @test_throws ArgumentError GadgetEditor.compute_payload(payload)
end
