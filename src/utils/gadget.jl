struct Gadget{T<:Real}
    ground_states::BitMatrix
    graph::SimpleGraph{Int}
    pins::Vector{Int}
    weights::Vector{T}
    pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}
end


function save_results_to_json(results::Vector{Gadget{T}}, file_path::String) where T
    # Convert each result to a serializable dictionary
    json_results = map(results) do res
        base_dict = Dict(
            "ground_states" => [vec(res.ground_states[i, :]) for i in 1:size(res.ground_states, 1)],
            "pins" => res.pins,
            "graph" => Dict(
                "nodes" => [Dict("id" => i, "weight" => res.weights[i]) for i in 1:length(res.weights)],
                "edges" => [Dict("source" => src(e), "target" => dst(e)) for e in Graphs.edges(res.graph)],
            )
        )

        # Add position data if it's a grid_gadget
        if !isnothing(res.pos)
            base_dict["graph"]["positions"] = [Dict("id" => i, "position" => [pos[1], pos[2]]) for (i, pos) in enumerate(res.pos)]
        end

        return base_dict
    end

    # Write to JSON file with pretty formatting
    open(file_path, "w") do io
        write(io, JSON3.write(json_results; pretty=true))
    end

    return file_path
end
