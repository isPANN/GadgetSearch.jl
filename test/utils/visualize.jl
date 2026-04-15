using Test
using GadgetSearch
using Graphs

function _make_sample_gadget()
    g = SimpleGraph(4)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 4)
    pos = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (1.0, 1.0)]
    weights = [1.0, 2.0, 3.0, 0.0]
    pins = [1, 2]
    tt = BitMatrix([1 0; 0 1])
    return GadgetSearch.Gadget(tt, g, pins, weights, pos)
end

@testset "plot_gadget: basic output formats" begin
    gadget = _make_sample_gadget()
    pdf1 = tempname() * ".pdf"
    svg1 = tempname() * ".svg"
    GadgetSearch.plot_gadget(gadget, pdf1)
    GadgetSearch.plot_gadget(gadget, svg1; show_weights=false)
    @test isfile(pdf1)
    @test isfile(svg1)
    rm(pdf1; force=true)
    rm(svg1; force=true)
end

@testset "plot_gadget: with options" begin
    gadget = _make_sample_gadget()
    pdf1 = tempname() * ".pdf"
    GadgetSearch.plot_gadget(gadget, pdf1;
        background_grid=true, shape="grid",
        show_weights=true, round_weights=true)
    @test isfile(pdf1)
    rm(pdf1; force=true)
end

@testset "plot_gadget: shape-aware grid" begin
    gadget = _make_sample_gadget()
    for s in ["grid", "KSG", "TLSG"]
        path = tempname() * ".svg"
        GadgetSearch.plot_gadget(gadget, path; shape=s, background_grid=true)
        @test isfile(path)
        rm(path; force=true)
    end
end

@testset "plot_graph: with and without positions" begin
    g = SimpleGraph(4)
    add_edge!(g, 1, 2); add_edge!(g, 2, 3); add_edge!(g, 3, 4)
    pos = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (1.0, 1.0)]
    pdf1 = tempname() * ".pdf"
    svg1 = tempname() * ".svg"
    GadgetSearch.plot_graph(g, pdf1; pos=nothing)
    GadgetSearch.plot_graph(g, svg1; pos=pos)
    @test isfile(pdf1)
    @test isfile(svg1)
    rm(pdf1; force=true)
    rm(svg1; force=true)
end

@testset "plot_single_gadget_mis" begin
    gadget = _make_sample_gadget()
    base_svg = tempname() * ".svg"
    GadgetSearch.plot_single_gadget_mis(gadget, base_svg)
    first_svg = replace(base_svg, ".svg" => "_1.svg")
    @test isfile(first_svg)
    for i in 1:5
        candidate = replace(base_svg, ".svg" => "_$(i).svg")
        isfile(candidate) && rm(candidate; force=true)
    end
end

@testset "_positions_to_layout" begin
    pts = GadgetSearch._positions_to_layout(
        [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)], 500, 20)
    @test length(pts) == 4
    # Collinear points
    pts2 = GadgetSearch._positions_to_layout([(0.0, 0.0), (1.0, 0.0)], 500, 20)
    @test length(pts2) == 2
end

@testset "_filter_zero_weights" begin
    gadget = _make_sample_gadget()
    graph, valid, old_to_new, weights, pins, pos, ew, el =
        GadgetSearch._filter_zero_weights(gadget)
    @test nv(graph) == 3  # vertex 4 (weight=0) removed
    @test length(weights) == 3
    @test all(w -> w != 0, weights)
end
