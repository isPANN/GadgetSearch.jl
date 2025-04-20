"""
    search_rules_parallel(graph_path, bit_num, gate_list, save_path;
                    pin_set::Vector{Int} = Int[],
                    # File processing parameters
                    max_file_size_mb::Int=30, split_size::Int=700_000,
                    # Search strategy parameters
                    start_idx::Int=0, end_idx::Int=0,
                    greedy::Bool=false, threshold::Int=0, max_samples::Int=0,
                    # Grid graph parameters
                    is_grid_graph::Bool=false, pos_data_path=nothing, grid_dims=(0,0),
                    # Optimizer Settings
                    optimizer=HiGHS.Optimizer, env=nothing,
                    # Parallel processing
                    num_threads::Int=Threads.nthreads())

Parallel version of search_rules that processes graph chunks in parallel.

# Arguments
- `graph_path`: Path to the graph database
- `bit_num`: Number of bits or vector of input/output bits
- `gate_list`: List of gate IDs to search for
- `save_path`: Path to save results

# Keyword Arguments
- `num_threads::Int=Threads.nthreads()`: Number of threads to use for parallel processing
- Other parameters are the same as search_rules

# Returns
- List of rule IDs for which no solution was found
"""
function search_rules_parallel(graph_path, bit_num, gate_list, save_path;
                    pin_set::Vector{Int} = Int[],
                    # File processing parameters
                    max_file_size_mb::Int=30, split_size::Int=700_000,
                    # Search strategy parameters
                    start_idx::Int=0, end_idx::Int=0,
                    greedy::Bool=false, threshold::Int=0, max_samples::Int=0,
                    # Grid graph parameters
                    is_grid_graph::Bool=false, pos_data_path=nothing, grid_dims=(0,0),
                    # Optimizer Settings
                    optimizer=HiGHS.Optimizer, env=nothing,
                    # Parallel processing parameters
                    batch_size::Int=50,
                    save_frequency::Int=10,  # 每找到多少个解决方案保存一次结果
                    log_frequency::Int=500   # 每处理多少个图打印一次日志
                    )

    # Determine graph type using traits
    graph_type = is_grid_graph ? GridGraph(pos_data_path, grid_dims) : GeneralGraph()

    # Create appropriate result vector based on graph type
    save_res = create_result_vector(graph_type)

    # Create search parameters object
    params = SearchParameters(
        pin_set=pin_set,
        max_file_size_mb=max_file_size_mb,
        split_size=split_size,
        start_idx=start_idx,
        end_idx=end_idx,
        greedy=greedy,
        threshold=threshold,
        max_samples=max_samples
    )

    # Prepare graph chunks for processing
    graph_chunks = if (filesize(graph_path) / (1024 * 1024)) > params.max_file_size_mb
        # Split the file into smaller parts and return as (index, path) pairs
        [(i, path) for (i, path) in enumerate(_split_large_file(graph_path, params.split_size))]
    else
        # Just use the original file with index 0
        [(0, graph_path)]
    end

    # 跟踪尚未找到解决方案的规则
    remaining_rules = Set(gate_list)

    # 确定搜索策略
    search_strategy = determine_search_strategy(bit_num)

    # 预先计算所有规则的基态
    rule_ground_states = Dict(rule_id => generic_rule(rule_id, bit_num) for rule_id in gate_list)

    # 创建线程安全的锁和计数器
    results_lock = ReentrantLock()
    rules_lock = ReentrantLock()

    # 创建线程本地结果存储，减少锁竞争
    thread_local_results = [Vector{AbstractGadget}() for _ in 1:Threads.nthreads()]
    thread_local_found_count = zeros(Int, Threads.nthreads())

    # 创建统计计数器
    processed_graphs = Threads.Atomic{Int}(0)
    found_solutions = Threads.Atomic{Int}(0)

    # 记录开始时间
    start_time = time()

    # 显示当前使用的线程数和其他参数
    @info "Starting parallel search with $(Threads.nthreads()) threads"
    @info "Batch size: $batch_size, Save frequency: $save_frequency, Log frequency: $log_frequency"

    # 预热函数，避免运行时编译引起的性能问题
    if !isempty(gate_list) && !isempty(graph_chunks)
        @info "Preheating functions..."
        dummy_graph = SimpleGraph(3)
        add_edge!(dummy_graph, 1, 2)
        dummy_mis, _ = find_maximal_independent_sets(dummy_graph)
        dummy_rule_id = first(gate_list)
        dummy_ground_states = rule_ground_states[dummy_rule_id]
        _check_grstates_in_candidate(dummy_mis, [1, 3], dummy_ground_states, false)
        @info "Preheating completed"
    end

    # 处理每个图块
    for (chunk_idx, chunk_path) in graph_chunks
        # 如果所有规则都已找到解决方案，提前退出
        isempty(remaining_rules) && break

        # 打印当前块信息
        chunk_start_time = time()
        @info "Processing chunk $(chunk_idx): $(chunk_path)"

        # Load graph dictionary once per chunk
        graph_dict = read_graph_dict(chunk_path)

        # 计算基础偏移量
        base_offset = chunk_idx == 0 ? 0 : (chunk_idx - 1) * params.split_size

        # 将图分成批次进行并行处理
        graph_batches = collect(enumerate(graph_dict))

        # 如果指定了起始和结束索引，过滤图
        if params.start_idx > 0 || params.end_idx > 0
            graph_batches = filter(pair -> begin
                i, _ = pair
                graph_id = base_offset + i
                (params.start_idx <= 0 || graph_id >= params.start_idx) &&
                (params.end_idx <= 0 || graph_id <= params.end_idx)
            end, graph_batches)
        end

        # 将图分成批次
        total_graphs = length(graph_batches)
        batch_indices = collect(1:batch_size:total_graphs)

        # 打印批次信息
        @info "Processing $(total_graphs) graphs in $(length(batch_indices)) batches"

        # 初始化每个线程的求解器环境，避免竞争
        thread_envs = [isnothing(env) ? nothing : deepcopy(env) for _ in 1:Threads.nthreads()]

        # 处理每个批次
        for batch_start in batch_indices
            # 如果所有规则都已找到解决方案，提前退出
            isempty(remaining_rules) && break

            # 计算批次结束索引
            batch_end = min(batch_start + batch_size - 1, total_graphs)
            batch = graph_batches[batch_start:batch_end]

            # 并行处理批次中的每个图
            Threads.@threads for (i, (gname, graph)) in batch
                # 当前线程ID
                thread_id = Threads.threadid()

                # 获取线程本地的求解器环境
                thread_env = thread_envs[thread_id]

                # 检查是否还有规则需要处理
                local_rules = lock(rules_lock) do
                    isempty(remaining_rules) ? Int[] : copy(remaining_rules)
                end

                # 如果没有规则需要处理，跳过当前图
                isempty(local_rules) && continue

                # 检查图是否连通，如果不连通则跳过
                Graphs.is_connected(graph) || continue

                # 计算原始图ID
                graph_id = _extract_numbers(gname) + base_offset

                # 更新处理的图计数
                Threads.atomic_add!(processed_graphs, 1)
                current_processed = processed_graphs[]

                # 打印进度信息
                if current_processed % log_frequency == 0
                    elapsed = time() - start_time
                    graphs_per_second = current_processed / elapsed
                    remaining_count = length(remaining_rules)
                    found_count = found_solutions[]
                    @info "Progress: processed $(current_processed) graphs in $(round(elapsed, digits=2))s ($(round(graphs_per_second, digits=2)) graphs/s), found $(found_count) solutions, $(remaining_count) rules remaining"
                end

                # 计算最大独立集（MIS）- 这是计算密集型操作，只需计算一次
                mis_result, mis_num = find_maximal_independent_sets(graph)

                # 根据搜索策略类型选择不同的处理方法
                if isa(search_strategy, LogicGateSearch)
                    input_num = bit_num[1]
                    output_num = bit_num[2]

                    # 如果MIS数量不足，跳过当前图
                    if mis_num < 2^(input_num)
                        continue
                    end

                    # 生成输入引脚候选
                    local pin_set = params.pin_set
                    if isempty(pin_set)
                        # 如果未提供引脚集，所有顶点都被视为潜在引脚
                        pin_set = collect(1:Graphs.nv(graph))
                        # 生成有效的输入引脚组合
                        input_candidates = _generate_pin_set(graph, mis_result, input_num)
                    else
                        # 如果提供了引脚集，则默认pin_set中的所有选择都有效
                        input_candidates = _generate_gate_input(pin_set, input_num)
                    end

                    # 如果没有有效的输入候选，跳过当前图
                    isempty(input_candidates) && continue

                    # 尝试每个输入候选
                    for candidate in input_candidates
                        # 检查是否还有规则需要处理
                        isempty(local_rules) && break

                        # 找到输出引脚的剩余元素
                        remain_elements = setdiff(pin_set, candidate)

                        # 尝试每种可能的输出引脚组合
                        for output_bits in permutations(remain_elements, output_num)
                            # 检查是否还有规则需要处理
                            isempty(local_rules) && break

                            # 组合输入和输出引脚
                            candidate_full = vcat(candidate, output_bits)

                            # 对每个剩余的规则进行处理
                            found_any = false

                            for rule_id in collect(local_rules)
                                # 检查规则是否已被其他线程处理
                                lock(rules_lock) do
                                    if !(rule_id in remaining_rules)
                                        delete!(local_rules, rule_id)
                                        return
                                    end
                                end

                                # 获取该规则的基态
                                ground_states = rule_ground_states[rule_id]

                                # 检查基态是否包含在MIS中
                                target_mis_indices_all = _check_grstates_in_candidate(mis_result, candidate_full, ground_states, params.greedy)
                                isempty(target_mis_indices_all) && continue

                                # 将target_mis_indices_all转换为_find_weight_new期望的格式
                                prefix_buckets = [collect(indices) for indices in target_mis_indices_all]

                                # 寻找权重 - 使用线程本地的求解器环境
                                weight = _find_weight_new(mis_result, prefix_buckets, optimizer, thread_env)
                                isempty(weight) && continue

                                # 如果找到解决方案
                                result = convert_to_gadget(graph_type, graph_id, graph, candidate_full, weight, ground_states, rule_id)

                                # 将结果存储在线程本地结果中，减少锁竞争
                                push!(thread_local_results[thread_id], result)
                                thread_local_found_count[thread_id] += 1
                                found_any = true

                                # 更新找到的解决方案计数
                                Threads.atomic_add!(found_solutions, 1)

                                # 从剩余规则中移除
                                lock(rules_lock) do
                                    delete!(remaining_rules, rule_id)
                                    delete!(local_rules, rule_id)
                                    println("Thread $(thread_id): Found solution for rule $(rule_id) on graph $(graph_id), $(length(remaining_rules)) rules remaining")
                                end

                                # 如果线程本地结果达到一定数量，就合并到全局结果并保存
                                if thread_local_found_count[thread_id] >= save_frequency
                                    lock(results_lock) do
                                        append!(save_res, thread_local_results[thread_id])
                                        save_results_to_json(save_res, save_path)
                                        # 清空线程本地结果
                                        empty!(thread_local_results[thread_id])
                                        thread_local_found_count[thread_id] = 0
                                    end
                                end
                            end

                            # 如果找到了解决方案，更新内存使用情况
                            if found_any && thread_id % 4 == 0  # 只有部分线程执行垃圾回收
                                GC.gc(false)  # 非阻塞垃圾回收，减少内存压力
                            end
                        end
                    end
                else  # 对于GenericConstraintSearch策略
                    # 生成所有引脚向量，即`pin_set`中`bit_num`个顶点的排列
                    local pin_set = params.pin_set
                    if isempty(pin_set)
                        pin_set = collect(1:Graphs.nv(graph))
                    end
                    @assert length(pin_set) >= bit_num
                    all_candidates = _generate_constraint_bit(pin_set, bit_num)

                    # 遍历所有可能的引脚向量
                    for candidate in all_candidates
                        # 检查是否还有规则需要处理
                        isempty(local_rules) && break

                        # 寴每个剩余的规则进行处理
                        found_any = false

                        for rule_id in collect(local_rules)
                            # 检查规则是否已被其他线程处理
                            lock(rules_lock) do
                                if !(rule_id in remaining_rules)
                                    delete!(local_rules, rule_id)
                                    return
                                end
                            end

                            # 获取该规则的基态
                            ground_states = rule_ground_states[rule_id]

                            # 检查在选择`candidate`的情况下，`ground_states`是否包含在MIS中
                            target_mis_indices_all = _check_grstates_in_candidate(mis_result, candidate, ground_states, params.greedy)
                            isempty(target_mis_indices_all) && continue

                            # 将target_mis_indices_all转换为_find_weight_new期望的格式
                            prefix_buckets = [collect(indices) for indices in target_mis_indices_all]

                            # 寻找权重 - 使用线程本地的求解器环境
                            weight = _find_weight_new(mis_result, prefix_buckets, optimizer, thread_env)
                            isempty(weight) && continue

                            # 如果找到解决方案
                            result = convert_to_gadget(graph_type, graph_id, graph, candidate, weight, ground_states, rule_id)

                            # 将结果存储在线程本地结果中，减少锁竞争
                            push!(thread_local_results[thread_id], result)
                            thread_local_found_count[thread_id] += 1
                            found_any = true

                            # 更新找到的解决方案计数
                            Threads.atomic_add!(found_solutions, 1)

                            # 从剩余规则中移除
                            lock(rules_lock) do
                                delete!(remaining_rules, rule_id)
                                delete!(local_rules, rule_id)
                                println("Thread $(thread_id): Found solution for rule $(rule_id) on graph $(graph_id), $(length(remaining_rules)) rules remaining")
                            end

                            # 如果线程本地结果达到一定数量，就合并到全局结果并保存
                            if thread_local_found_count[thread_id] >= save_frequency
                                lock(results_lock) do
                                    append!(save_res, thread_local_results[thread_id])
                                    save_results_to_json(save_res, save_path)
                                    # 清空线程本地结果
                                    empty!(thread_local_results[thread_id])
                                    thread_local_found_count[thread_id] = 0
                                end
                            end
                        end

                        # 如果找到了解决方案，更新内存使用情况
                        if found_any && thread_id % 4 == 0  # 只有部分线程执行垃圾回收
                            GC.gc(false)  # 非阻塞垃圾回收，减少内存压力
                        end
                    end
                end
            end
        end
    end

    # 合并所有线程本地结果到全局结果并保存
    lock(results_lock) do
        for thread_id in 1:Threads.nthreads()
            if !isempty(thread_local_results[thread_id])
                append!(save_res, thread_local_results[thread_id])
                empty!(thread_local_results[thread_id])
            end
        end

        # 保存最终结果
        if !isempty(save_res)
            save_results_to_json(save_res, save_path)
        end
    end

    # 打印最终统计信息
    elapsed = time() - start_time
    total_processed = processed_graphs[]
    total_found = found_solutions[]

    @info "Search completed in $(round(elapsed, digits=2)) seconds"
    @info "Processed $(total_processed) graphs ($(round(total_processed/elapsed, digits=2)) graphs/s)"
    @info "Found $(total_found) solutions for $(length(gate_list) - length(remaining_rules)) rules"
    @info "$(length(remaining_rules)) rules remain unsolved"

    # 返回未找到解决方案的规则
    return collect(remaining_rules)
end
