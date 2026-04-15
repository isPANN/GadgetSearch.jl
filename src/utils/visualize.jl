# ----- Internal helpers -----

function _with_drawing(drawfn::Function, save_path::AbstractString, w::Integer, h::Integer)
    Drawing(w, h, save_path)
    try
        origin()
        drawfn()
    finally
        finish()
    end
    return save_path
end

function _suffix_path(path::AbstractString, suffix::AbstractString)
    base, ext = splitext(path)
    return string(base, suffix, ext)
end

function _positions_to_layout(pos::Vector{Tuple{Float64, Float64}}, plot_size::Int, margin::Int; preserve_aspect_ratio::Bool=true)
    x_vals = [p[1] for p in pos]
    y_vals = [p[2] for p in pos]
    x_center = (minimum(x_vals) + maximum(x_vals)) / 2
    y_center = (minimum(y_vals) + maximum(y_vals)) / 2
    xc = x_vals .- x_center
    yc = y_vals .- y_center
    allowed = plot_size / 2 - margin
    x_ext = max(maximum(abs.(xc)), eps(Float64))
    y_ext = max(maximum(abs.(yc)), eps(Float64))
    sx = allowed / x_ext
    sy = allowed / y_ext
    if preserve_aspect_ratio
        s = min(sx, sy)
        sx, sy = s, s
    end
    return [Point(x * sx, -y * sy) for (x, y) in zip(xc, yc)]
end

function _generate_vertex_color(weights::AbstractVector, discrete_color_scheme, continuous_color_scheme)
    if any(w < 1 || w > length(discrete_color_scheme) || typeof(w) != Int for w in weights)
        min_w, max_w = minimum(weights), maximum(weights)
        cmap = get(continuous_color_scheme, range(0, 1, length=length(weights)))
        if max_w == min_w
            mid = clamp(round(Int, length(cmap)/2), 1, length(cmap))
            return fill(cmap[mid], length(weights))
        end
        return [cmap[clamp(floor(Int, (w - min_w) / (max_w - min_w) * (length(cmap) - 1) + 1), 1, length(cmap))] for w in weights]
    else
        return [discrete_color_scheme[w] for w in weights]
    end
end

function _filter_zero_weights(gadget::Gadget)
    valid = findall(w -> w != 0, gadget.weights)
    if length(valid) == nv(gadget.graph)
        return (gadget.graph, valid, Dict(v => v for v in valid),
                gadget.weights, gadget.pins, gadget.pos,
                gadget.edge_weights, gadget.edge_list)
    end
    old_to_new = Dict(v => i for (i, v) in enumerate(valid))
    new_graph = SimpleGraph(length(valid))
    has_ew = !isempty(gadget.edge_weights)
    new_edge_weights = Float64[]
    new_edge_list = Tuple{Int,Int}[]
    if has_ew
        for (idx, (s, d)) in enumerate(gadget.edge_list)
            if s in valid && d in valid
                ns, nd = old_to_new[s], old_to_new[d]
                add_edge!(new_graph, ns, nd)
                push!(new_edge_list, (ns, nd))
                push!(new_edge_weights, gadget.edge_weights[idx])
            end
        end
    else
        for e in edges(gadget.graph)
            s, d = Graphs.src(e), Graphs.dst(e)
            if s in valid && d in valid
                add_edge!(new_graph, old_to_new[s], old_to_new[d])
            end
        end
    end
    new_weights = gadget.weights[valid]
    new_pins = [old_to_new[p] for p in gadget.pins if p in valid]
    new_pos = gadget.pos !== nothing ? gadget.pos[valid] : nothing
    return (new_graph, valid, old_to_new,
            new_weights, new_pins, new_pos,
            new_edge_weights, new_edge_list)
end

# ----- Public API -----

function plot_gadget(gadget::Gadget, save_path::String;
                     plot_size=400, margin=30,
                     preserve_aspect_ratio=true,
                     show_weights=true,
                     show_edge_weights=true,
                     round_weights=false,
                     discrete_color_scheme=ColorSchemes.seaborn_bright,
                     continuous_color_scheme=ColorSchemes.viridis)
    graph, _, _, weights, pins, pos, edge_weights, edge_list = _filter_zero_weights(gadget)
    pos === nothing && error("Gadget has no positions; cannot plot")
    pts = _positions_to_layout(pos, plot_size, margin; preserve_aspect_ratio=preserve_aspect_ratio)
    node_colors = _generate_vertex_color(weights, discrete_color_scheme, continuous_color_scheme)
    has_ew = !isempty(edge_weights)
    _with_drawing(save_path, plot_size, plot_size) do
        background("white")
        sethue("black")
        # Edge weight map for QUBO
        edge_weight_map = Dict{Tuple{Int,Int}, Float64}()
        if has_ew
            for (idx, (s, d)) in enumerate(edge_list)
                key = s < d ? (s, d) : (d, s)
                edge_weight_map[key] = edge_weights[idx]
            end
        end
        # Draw edges with labels if QUBO
        if has_ew && show_edge_weights
            for e in edges(graph)
                s, d = Graphs.src(e), Graphs.dst(e)
                key = s < d ? (s, d) : (d, s)
                sethue("black"); setline(2)
                line(pts[s], pts[d], :stroke)
                if haskey(edge_weight_map, key)
                    w = edge_weight_map[key]
                    mid = Point((pts[s].x + pts[d].x) / 2, (pts[s].y + pts[d].y) / 2)
                    dx = pts[d].x - pts[s].x; dy = pts[d].y - pts[s].y
                    len = sqrt(dx^2 + dy^2)
                    perp = Point(-dy / len, dx / len)
                    lp = mid + perp * 12
                    wt = round_weights ? string(round(Int, w)) : string(w)
                    sethue("white"); fontsize(11)
                    tw = textextents(wt)[3]; th = textextents(wt)[4]
                    box(lp, tw + 4, th + 4, :fill)
                    sethue("blue"); fontsize(11)
                    text(wt, lp, halign=:center, valign=:middle)
                end
            end
        end
        drawgraph(graph,
                  vertexfunction=(v, c) -> begin
                      pin_idx = findfirst(==(v), pins)
                      sethue(node_colors[v])
                      translate(pts[v])
                      circle(O, 10, :fill)
                      if show_weights
                          sethue("white"); fontsize(15)
                          wt = round_weights ? "$(round(Int, weights[v]))" : "$(weights[v])"
                          text(wt, O, halign=:center, valign=:middle)
                      end
                      if pin_idx !== nothing
                          sethue("red")
                          dy = pts[v].y >= 0 ? 20 : -20
                          text("PIN $pin_idx", Point(0, dy), halign=:center, valign=:middle)
                      end
                  end,
                  edgestrokeweights=has_ew && show_edge_weights ? 0 : 2,
                  layout=pts)
    end
    println("Drawing saved as $save_path")
end

function plot_graph(g::SimpleGraph, save_path::String;
                    pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}=nothing,
                    plot_size=400, margin=30,
                    vertex_size=10, vertex_label_size=15, edge_width=2)
    _with_drawing(save_path, plot_size, plot_size) do
        background("white")
        layout = pos === nothing ? stress : _positions_to_layout(pos, plot_size, margin)
        sethue("black")
        drawgraph(g,
                  vertexlabels=collect(1:nv(g)),
                  vertexshapesizes=vertex_size,
                  vertexlabelfontsizes=vertex_label_size,
                  edgestrokeweights=edge_width,
                  layout=layout)
    end
    println("Drawing saved as $save_path")
end

function plot_single_gadget_mis(gadget::Gadget, save_path::String;
                                plot_size=400, margin=30)
    masks, n = find_maximal_independent_sets(gadget.graph)
    for i in 1:n
        save_single_path = _suffix_path(save_path, "_" * string(i))
        _with_drawing(save_single_path, plot_size, plot_size) do
            m = masks[i]
            colors = [(((m >> (v - 1)) & 0x1) == 1) ? RGB(1, 0, 0) : RGB(0, 0, 0) for v in 1:nv(gadget.graph)]
            drawgraph(gadget.graph,
                      vertexlabels=collect(1:nv(gadget.graph)),
                      vertexshapesizes=18,
                      vertexlabelfontsizes=22,
                      vertexfillcolors=colors,
                      edgestrokeweights=3,
                      layout=spring)
        end
        println("Drawing saved as $save_single_path")
    end
end
