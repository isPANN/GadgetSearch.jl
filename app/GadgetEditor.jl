module GadgetEditor

using GadgetSearch
using Graphs
using HTTP
using JSON3

const APP_DIR = @__DIR__

function gadget_from_payload(payload)
    nodes = payload.nodes
    isempty(nodes) && throw(ArgumentError("Place at least one vertex on the lattice"))
    weight_mode = String(payload.weight_mode)
    weight_mode in ("weighted", "unweighted") ||
        throw(ArgumentError("weight_mode must be weighted or unweighted"))
    if weight_mode == "weighted" && length(nodes) > GadgetSearch.MAX_SUPPORTED_VERTICES
        throw(ArgumentError(
            "Weighted ground-state enumeration supports at most " *
            "$(GadgetSearch.MAX_SUPPORTED_VERTICES) vertices",
        ))
    end

    ids = String[String(node.id) for node in nodes]
    length(unique(ids)) == length(ids) || throw(ArgumentError("Vertex IDs must be unique"))
    id_to_index = Dict(id => index for (index, id) in enumerate(ids))

    lattice_coordinates = Tuple{Int, Int}[(Int(node.q), Int(node.r)) for node in nodes]
    positions = GadgetSearch.get_physical_positions(Triangular(), lattice_coordinates)
    weights = weight_mode == "weighted" ?
        Float64[Float64(node.weight) for node in nodes] :
        ones(Float64, length(nodes))
    weight_mode == "weighted" && !all(>(0), weights) &&
        throw(ArgumentError("Weighted mode requires positive vertex weights"))

    graph = GadgetSearch.unit_disk_graph(positions, get_radius(Triangular()))

    pin_ids = String[String(id) for id in payload.pins]
    all(id -> haskey(id_to_index, id), pin_ids) ||
        throw(ArgumentError("Pins must reference vertices on the lattice"))
    length(unique(pin_ids)) == length(pin_ids) ||
        throw(ArgumentError("A vertex cannot be used as a pin more than once"))
    pins = Int[id_to_index[id] for id in pin_ids]

    constraint = TruthTableConstraint(BitMatrix(falses(1, length(pins))))
    return (
        Gadget(RydbergModel, constraint, graph, pins, weights, positions),
        ids,
        weight_mode,
    )
end

function compute_payload(payload)
    String(payload.model) == "rydberg" ||
        throw(ArgumentError("The editor currently supports the Rydberg / MIS model"))

    gadget, ids, weight_mode = gadget_from_payload(payload)
    edge_data = [
        (source=ids[src(edge)], target=ids[dst(edge)]) for edge in edges(gadget.graph)
    ]

    if weight_mode == "unweighted"
        boundary_count = length(gadget.pins)
        tensor_result = calculate_reduced_alpha_tensor(gadget.graph, gadget.pins)
        tensor = boundary_count == 0 ? [tensor_result] : vec(tensor_result)
        entries = [
            (
                configuration=join(
                    Int[((index - 1) >> (bit - 1)) & 1 for bit in 1:boundary_count]
                ),
                value=isfinite(value) ? value : "-Inf",
            )
            for (index, value) in enumerate(tensor)
        ]
        return (
            operation="reduced_alpha_tensor",
            model="Unweighted MIS",
            weight_mode=weight_mode,
            vertex_count=nv(gadget.graph),
            edge_count=ne(gadget.graph),
            boundary_count=boundary_count,
            edges=edge_data,
            tensor=entries,
        )
    end

    report = analyze_gadget(gadget; model=RydbergModel)
    observed = [join(state.pins) for state in report.ground_states]
    ground_states = [
        (
            state_index=state.state_index,
            pins=join(state.pins),
            configuration=join(state.configuration),
            occupied=ids[findall(==(1), state.configuration)],
        )
        for state in report.ground_states
    ]

    return (
        operation="ground_states",
        model=report.model,
        weight_mode=weight_mode,
        vertex_count=nv(gadget.graph),
        edge_count=ne(gadget.graph),
        state_count=report.state_count,
        max_energy=report.max_energy,
        observed=observed,
        edges=edge_data,
        ground_states=ground_states,
    )
end

json_response(status::Int, body) = HTTP.Response(
    status,
    ["Content-Type" => "application/json; charset=utf-8"],
    JSON3.write(body),
)

function static_response(filename::String, content_type::String)
    return HTTP.Response(
        200,
        ["Content-Type" => content_type],
        read(joinpath(APP_DIR, filename)),
    )
end

function handler(request::HTTP.Request)
    path = HTTP.URI(request.target).path
    if request.method == "POST" && path == "/api/compute"
        try
            return json_response(200, compute_payload(JSON3.read(request.body)))
        catch error
            return json_response(422, (error=sprint(showerror, error),))
        end
    elseif request.method == "GET" && (path == "/" || path == "/index.html")
        return static_response("index.html", "text/html; charset=utf-8")
    elseif request.method == "GET" && path == "/styles.css"
        return static_response("styles.css", "text/css; charset=utf-8")
    elseif request.method == "GET" && path == "/app.js"
        return static_response("app.js", "text/javascript; charset=utf-8")
    end
    return json_response(404, (error="Not found",))
end

function serve(; host="127.0.0.1", port=8080)
    println("Gadget Editor is running at http://$(host):$(port)")
    HTTP.serve(handler, host, port)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    GadgetEditor.serve()
end
