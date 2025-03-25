# TODO to be refactored
function search_single_constraint_udg(graph_path::String, bit_num::Vector{Int}, ground_states::Vector{Int}; gate_id::Int = 0, max_file_size::Int = 1_000_000)::Union{GadgetSolution, Nothing}
    if filesize(graph_path) > max_file_size
        _process_large_file(graph_path, bit_num, ground_states, gate_id; pin_set=[i for i in 1:sum(bit_num)])
    else
        _execute_graph_search(graph_path, bit_num, ground_states, gate_id; pin_set=[i for i in 1:sum(bit_num)])
    end
end

function search_single_constraint_udg(graph_path::String, bit_num::Vector{Int}, gate_id::Int)
    @assert length(bit_num) == 2
    input_num = bit_num[1]; output_num = bit_num[2];
    return search_single_constraint_udg(graph_path, bit_num, generic_gate(gate_id, input_num, output_num); gate_id=gate_id)
end

function search_gates_udg(graph_path, bit_num, gate_list, save_path)
    non_results = []
    save_res = []
    for gate_id in gate_list
        res = search_single_constraint_udg(graph_path, bit_num, gate_id)
        if !isnothing(res)
            push!(save_res, res)
            save_results_to_json(save_res, save_path)
            res = nothing
        else
            push!(non_results, gate_id)
        end
    end
    return non_results
end