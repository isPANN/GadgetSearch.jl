using Test
using GadgetSearch
using Graphs

"""
    _make_sample_gadget()

Create a minimal Gadget suitable for visualization tests.
"""
function _make_sample_gadget()
    g = SimpleGraph(4)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 4)
    pos = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (1.0, 1.0)]
    weights = [1.0, 2.0, 3.0, 0.0]  # include a zero to exercise filtering in plot_single_gadget_new
    pins = [1, 2]
    tt = BitMatrix([1 0; 0 1])
    return GadgetSearch.Gadget(tt, g, pins, weights, pos)
end

@testset "visualize: output formats and basics" begin
    gadget = _make_sample_gadget()
    g = gadget.graph

    # plot_single_gadget to PDF and SVG
    pdf1 = tempname() * ".pdf"
    svg1 = tempname() * ".svg"
    GadgetSearch.plot_single_gadget(gadget, pdf1)
    GadgetSearch.plot_single_gadget(gadget, svg1)
    @test isfile(pdf1)
    @test isfile(svg1)
    rm(pdf1; force=true)
    rm(svg1; force=true)

    # plot_single_gadget_new with/without options
    pdf2 = tempname() * ".pdf"
    svg2 = tempname() * ".svg"
    GadgetSearch.plot_gadget(gadget, pdf2; background_grid=true, show_weights=true, round_weights=true)
    GadgetSearch.plot_gadget(gadget, svg2; background_grid=false, show_weights=false)
    @test isfile(pdf2)
    @test isfile(svg2)
    rm(pdf2; force=true)
    rm(svg2; force=true)

    # plot_graph with and without positions
    pdf3 = tempname() * ".pdf"
    svg3 = tempname() * ".svg"
    GadgetSearch.plot_graph(g, pdf3; pos=nothing)
    GadgetSearch.plot_graph(g, svg3; pos=gadget.pos)
    @test isfile(pdf3)
    @test isfile(svg3)
    rm(pdf3; force=true)
    rm(svg3; force=true)

    # plot_single_gadget_mis: ensure it writes multiple files with suffixes
    base_svg = tempname() * ".svg"
    GadgetSearch.plot_single_gadget_mis(gadget, base_svg)
    # Expect at least one file with _1 suffix
    first_svg = replace(base_svg, ".svg" => "_1.svg")
    @test isfile(first_svg)
    # Clean up a few possible outputs
    for i in 1:5
        candidate = replace(base_svg, ".svg" => "_$(i).svg")
        isfile(candidate) && rm(candidate; force=true)
    end
end
