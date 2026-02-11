# α-Tensor 实现指南（给 AI 的详细说明）

## 🎯 任务目标

在 **GadgetSearch.jl** 中实现 **α-tensor 验证模式**，基于论文：
> Liu et al., "Computer-assisted gadget design and problem reduction of unweighted maximum independent set", PRX Quantum 4, 010316 (2023)

---

## ❓ 核心设计问题与答案

### 问题 1：TruthTableConstraint 和 α-tensor 的关系

**答案**：它们是互补的概念，但可以相互转换。

#### 概念对应

```
TruthTableConstraint (真值表)
  ↓ 定义
【pins 上的目标配置】
  例如：OR gate 的真值表定义了 pins 上哪些配置应该是基态
  
α-tensor (约简张量)
  ↓ 编码
【整个图的 MIS 结构】
  例如：对于每个 pin 配置，内部最大独立集的大小
```

#### 转换关系

```julia
# 1. 从真值表 → 目标状态（MIS bitmask）
truth_table = BitMatrix([
    0 0 0;  # pin1=0, pin2=0 → output=0
    1 0 1;  # pin1=1, pin2=0 → output=1
    ...
])
# 对应的基态 MIS（例如）：{pin1, output}, {pin2, output}, ...
target_states = [0b101, 0b110, ...]  # 二进制表示

# 2. 从目标状态 → pattern 的 α-tensor
function infer_pattern_alpha(target_states, pins)
    alpha = Dict{UInt32, Int}()
    for state in target_states
        pin_config = extract_pin_config(state, pins)
        interior_size = count_ones(state) - count_ones(pin_config)
        alpha[pin_config] = max(get(alpha, pin_config, 0), interior_size)
    end
    return alpha
end

# 3. gadget 的 α-tensor（直接从图计算）
alpha_gadget = compute_reduced_alpha_tensor(graph, pins)

# 4. 验证等价性
is_valid = (alpha_gadget == alpha_pattern + constant)
```

#### 关键理解

- **TruthTableConstraint**：用户输入，定义逻辑功能
- **Target states**：从真值表推导出的具体 MIS
- **Pattern α-tensor**：从 target states 推断出的抽象特征
- **Gadget α-tensor**：从图结构计算出的实际特征
- **验证**：两个 α-tensor 是否相差一个常数

---

### 问题 2：α-tensor 的实现方式

**答案**：选择 **(A) 直接枚举独立集（暴力法）**

#### 为什么选择暴力法？

1. **GadgetSearch.jl 已有基础设施**
   ```julia
   # 已经实现并优化过
   find_maximal_independent_sets(graph)  # 返回所有 MIS
   ```

2. **问题规模小**
   - 最多 32 个顶点（UInt32 bitmask 限制）
   - Pin 数通常 ≤ 4（最多 16 个配置）
   - 暴力法完全可行

3. **tropical tensor network 的复杂性**
   - 需要实现 max-plus semiring
   - 需要张量收缩算法
   - 对小规模问题收益不大

#### 具体算法（我的实现）

```julia
function compute_reduced_alpha_tensor(graph, pins)
    n_pins = length(pins)
    interior = setdiff(1:nv(graph), pins)
    alpha = Dict{UInt32, Int}()
    
    # 遍历所有 2^n_pins 个边界配置
    for boundary_mask in 0:(2^n_pins - 1)
        # 1. 确定哪些 pins 被选中
        selected_pins = [pins[i] for i in 1:n_pins if (boundary_mask >> (i-1)) & 1 == 1]
        
        # 2. 找出被选中 pins 的邻居（内部顶点）
        forbidden = Set()
        for pin in selected_pins
            for neighbor in neighbors(graph, pin)
                if neighbor ∈ interior
                    push!(forbidden, neighbor)
                end
            end
        end
        
        # 3. 在剩余的内部顶点上找最大独立集
        feasible_interior = setdiff(interior, forbidden)
        if isempty(feasible_interior)
            alpha[UInt32(boundary_mask)] = 0
        else
            subgraph, _ = induced_subgraph(graph, feasible_interior)
            mis_masks, _ = find_maximal_independent_sets(subgraph)
            alpha[UInt32(boundary_mask)] = maximum(count_ones, mis_masks; init=0)
        end
    end
    
    return alpha
end
```

#### mis_compactify! 的实现

**我没有实现 mis_compactify!**，因为不需要。原因：

- `mis_compactify!` 用于压缩 α-tensor（去除被支配的配置）
- 在验证阶段，我们只需要检查等价性，不需要压缩
- 如果要实现：

```julia
# 伪代码
function mis_compactify!(alpha::Dict{UInt32, Int})
    # 遍历所有配置对
    for (config1, val1) in alpha
        for (config2, val2) in alpha
            if config1 != config2
                # 如果 config1 ⊆ config2（子集关系）
                if config1 & config2 == config1
                    # 且 val1 <= val2（被支配）
                    if val1 <= val2
                        delete!(alpha, config1)  # 删除被支配的
                    end
                end
            end
        end
    end
end
```

但这不是必需的，因为等价性检查不受压缩影响。

---

### 问题 3：目标 pattern 的定义

**答案**：Pattern 是从 TruthTableConstraint 动态推断的，不需要预定义。

#### 设计选择

**我的方案**：动态推断 pattern

```julia
# 用户输入：真值表
constraint = TruthTableConstraint(BitMatrix([...]))

# 搜索时：
# 1. 从真值表生成目标 MIS 状态
target_states = generate_target_states(constraint, pins)

# 2. 从目标状态推断 pattern 的 α-tensor
alpha_pattern = infer_pattern_alpha(target_states, pins)

# 3. 对每个候选 gadget：
#    - 计算 alpha_gadget
#    - 检查 alpha_gadget 与 alpha_pattern 是否等价
```

#### 与预定义 pattern 的对比

| 方案 | 优点 | 缺点 |
|------|------|------|
| **动态推断**（我的实现） | ✅ 灵活，支持任意真值表<br>✅ 不需要维护 pattern 库<br>✅ 与现有 API 一致 | ⚠️ 可能找到非标准结构 |
| **预定义 pattern** | ✅ 保证找到的是标准 gadget<br>✅ 可以复用论文中的 pattern | ⚠️ 需要手动定义所有 pattern<br>⚠️ 限制了搜索空间 |

#### 实现细节

```julia
function infer_pattern_alpha(target_states::Vector{UInt32}, pins::Vector{Int})
    alpha = Dict{UInt32, Int}()
    
    for state in target_states
        # 提取 pin 配置
        pin_config = UInt32(0)
        for (i, pin) in enumerate(pins)
            if (state >> (pin - 1)) & 0x1 == 1
                pin_config |= UInt32(1) << (i - 1)
            end
        end
        
        # 计算内部顶点数
        total_size = count_ones(state)
        n_selected_pins = count_ones(pin_config)
        interior_size = total_size - n_selected_pins
        
        # 更新 alpha（取最大值）
        if haskey(alpha, pin_config)
            alpha[pin_config] = max(alpha[pin_config], interior_size)
        else
            alpha[pin_config] = interior_size
        end
    end
    
    return alpha
end
```

#### 关键点

- **不需要**用户指定 pattern 类型（CROSS、WIRE 等）
- TruthTableConstraint 已经隐式定义了 pattern
- α-tensor 是 pattern 的数学表征，自动推断

---

### 问题 4：King's subgraph 嵌入检查

**答案**：选择 **(A) 直接用现有的 UDG 生成器限制搜索空间**

#### 实现策略

```julia
# 1. 数据集生成阶段：只生成 King's subgraph
generate_full_grid_udg(Triangular(), 3, 3; path="dataset.g6")
# ↑ 这已经保证了所有图都是 King's lattice 的子图

# 2. 搜索阶段：在这个受限数据集中搜索
loader = GraphLoader("dataset.g6")  # 只包含有效的 King's subgraph

# 3. 验证阶段：用 α-tensor 检查
result = verify_gadget_via_alpha_tensor(graph, pins, target_states)
```

#### 为什么不做运行时检查？

| 方案 | 优点 | 缺点 | 是否采用 |
|------|------|------|----------|
| (A) 限制搜索空间 | ✅ 简单高效<br>✅ 利用现有代码 | ⚠️ 需要预生成数据集 | ✅ **采用** |
| (B) 同构检查 | ✅ 支持任意图输入 | ❌ NP-complete 问题<br>❌ 实现复杂 | ❌ 不采用 |
| (C) 优化 loss function | ✅ 理论完整 | ❌ 需要优化器<br>❌ 违背"无优化器"目标 | ❌ 不采用 |

#### King's lattice 的特性

```julia
# Triangular lattice (论文使用的)
# 每个顶点最多 6 个邻居
# 可以嵌入到平面上，保持单位圆盘性质

# Square lattice with diagonals (King's graph)
# 每个顶点最多 8 个邻居
# 与 Triangular 质量因子 Q 不同
```

GadgetSearch.jl 的 `generate_full_grid_udg` 支持两种：

```julia
generate_full_grid_udg(Triangular(), 3, 3)  # 三角晶格
generate_full_grid_udg(Square(), 3, 3)      # 方形晶格（King's graph）
```

---

### 问题 5：架构集成方式

**答案**：选择 **(C) 两者并存，用 dispatch 区分**

#### 架构设计

```julia
# 1. 定义新的 EnergyModel
abstract type EnergyModel end
struct RydbergModel <: EnergyModel end
struct QUBOModel <: EnergyModel end
struct RydbergUnweightedModel <: EnergyModel end  # 简单基数检查
struct AlphaTensorMode <: EnergyModel end         # α-tensor 验证

# 2. 为每个模式实现 _find_weights
function _find_weights(::Type{RydbergModel}, ...)
    # 使用整数规划找权重
    # 需要 optimizer
end

function _find_weights(::Type{RydbergUnweightedModel}, ...)
    # 简单的基数检查（所有权重=1）
    # 不需要 optimizer
    target_energy = count_ones(target_states[1])
    for s in target_states
        count_ones(s) == target_energy || return nothing
    end
    for s in wrong_states
        count_ones(s) < target_energy || return nothing
    end
    return ones(Float64, vertex_num)
end

function _find_weights(::Type{AlphaTensorMode}, ...)
    # α-tensor 验证
    # 不需要 optimizer
    result = verify_gadget_via_alpha_tensor(graph, pins, target_states)
    return result === nothing ? nothing : result[1]  # 返回权重
end

# 3. 统一的搜索接口
function search_gadgets(
    ::Type{M},
    loader::GraphLoader,
    constraints::Vector{C};
    optimizer=nothing,  # 对 AlphaTensorMode 可选
    ...
) where {M <: EnergyModel, C <: GadgetConstraint}
    # optimizer 检查
    if optimizer === nothing && !(M <: Union{RydbergUnweightedModel, AlphaTensorMode})
        error("Optimizer required for $(M)")
    end
    
    # 统一的搜索流程
    for graph in loader
        result = _find_weights(M, ...)  # dispatch 到具体实现
        if result !== nothing
            push!(results, Gadget(...))
        end
    end
end
```

#### 为什么选择 dispatch？

| 方案 | 优点 | 缺点 |
|------|------|------|
| (A) 替换 _find_weights | ❌ 破坏现有功能 | ❌ 不兼容 |
| (B) 新建框架 | ✅ 完全独立 | ❌ 代码重复<br>❌ 用户困惑 |
| (C) Dispatch | ✅ 代码复用<br>✅ 统一接口<br>✅ 易于扩展 | ⚠️ 需要理解 Julia dispatch |

#### 关键实现点

```julia
# src/core/alpha_tensor.jl
# 独立的 α-tensor 函数，可以单独使用
compute_reduced_alpha_tensor(graph, pins)
check_alpha_equivalence(α1, α2)
verify_gadget_via_alpha_tensor(graph, pins, target_states)

# src/core/search.jl  
# 集成到搜索框架
struct AlphaTensorMode <: EnergyModel end

function _find_weights(::Type{AlphaTensorMode}, ...)
    # 调用 alpha_tensor.jl 的函数
    verify_gadget_via_alpha_tensor(...)
end
```

---

## 📋 完整实现步骤

### Step 1: 创建 alpha_tensor.jl

```julia
# src/core/alpha_tensor.jl

"""
计算约简 α-tensor
"""
function compute_reduced_alpha_tensor(
    graph::SimpleGraph{Int},
    pins::Vector{Int}
)
    n_pins = length(pins)
    interior = setdiff(1:nv(graph), pins)
    alpha = Dict{UInt32, Int}()
    
    # 遍历所有边界配置
    for boundary_mask in 0:(2^n_pins - 1)
        # ... (见问题2的代码)
        alpha[UInt32(boundary_mask)] = max_interior_mis_size
    end
    
    return alpha
end

"""
检查两个 α-tensor 是否等价
"""
function check_alpha_equivalence(
    α1::Dict{UInt32, Int},
    α2::Dict{UInt32, Int}
)
    # 键必须相同
    Set(keys(α1)) == Set(keys(α2)) || return (false, nothing)
    
    # 计算常数差
    first_config = first(keys(α1))
    c = α2[first_config] - α1[first_config]
    
    # 验证所有配置的差都是 c
    for config in keys(α1)
        if α2[config] - α1[config] != c
            return (false, nothing)
        end
    end
    
    return (true, c)
end

"""
从目标状态推断 pattern 的 α-tensor
"""
function infer_pattern_alpha(
    target_states::Vector{UInt32},
    pins::Vector{Int}
)
    # ... (见问题3的代码)
end

"""
完整的 gadget 验证流程
"""
function verify_gadget_via_alpha_tensor(
    graph::SimpleGraph{Int},
    pins::Vector{Int},
    target_states::Vector{UInt32}
)
    # 1. 推断 pattern 的 α-tensor
    α_pattern = infer_pattern_alpha(target_states, pins)
    
    # 2. 计算 gadget 的 α-tensor
    α_gadget = compute_reduced_alpha_tensor(graph, pins)
    
    # 3. 检查等价性
    is_equiv, overhead = check_alpha_equivalence(α_pattern, α_gadget)
    
    if !is_equiv
        return nothing
    end
    
    # 返回均匀权重和开销
    weights = ones(Float64, nv(graph))
    return (weights, overhead)
end
```

### Step 2: 修改 search.jl

```julia
# src/core/search.jl

# 添加新类型
struct AlphaTensorMode <: EnergyModel end

# 添加 dispatch
function get_state_space(::Type{AlphaTensorMode}, graph::SimpleGraph{Int})
    return find_maximal_independent_sets(graph)  # MIS state space
end

function _find_weights(
    ::Type{AlphaTensorMode},
    vertex_num::Int,
    edge_list::Vector{Tuple{Int,Int}},
    pin_set::Vector{Int},
    target_states::Vector{UInt32},
    wrong_states::Vector{UInt32},
    optimizer,  # 不使用
    env,
    objective,
    allow_defect::Bool,
    graph::SimpleGraph{Int},
    check_connectivity::Bool=true
)
    result = verify_gadget_via_alpha_tensor(graph, pin_set, target_states)
    
    if result === nothing
        return nothing
    end
    
    weights, overhead = result
    @info "found a valid AlphaTensor solution (overhead = $overhead)"
    return weights
end

# 修改 search_gadgets：optimizer 可选
function search_gadgets(
    ::Type{M},
    loader::GraphLoader,
    constraints::Vector{C};
    optimizer=nothing,  # 改为可选
    ...
) where {M <: EnergyModel, C <: GadgetConstraint}
    # 只有需要的模式才检查 optimizer
    if optimizer === nothing && !(M <: Union{RydbergUnweightedModel, AlphaTensorMode})
        error("Optimizer required for $(M)")
    end
    
    # ... 其余代码不变
end
```

### Step 3: 导出新功能

```julia
# src/GadgetSearch.jl

include("core/alpha_tensor.jl")  # 添加这行

export AlphaTensorMode  # 导出类型
export compute_reduced_alpha_tensor  # 导出函数
export check_alpha_equivalence
export verify_gadget_via_alpha_tensor
```

### Step 4: 编写测试

```julia
# test/core/alpha_tensor_test.jl

using Test, GadgetSearch, Graphs

@testset "α-Tensor Computation" begin
    # 测试简单图
    g = path_graph(3)
    pins = [2]
    α = compute_reduced_alpha_tensor(g, pins)
    
    @test α[UInt32(0)] == 2  # pin 2 不选 → {1,3}
    @test α[UInt32(1)] == 0  # pin 2 选中 → 邻居都被阻止
end

@testset "Equivalence Check" begin
    α1 = Dict(UInt32(0) => 2, UInt32(1) => 0)
    α2 = Dict(UInt32(0) => 5, UInt32(1) => 3)
    
    is_equiv, c = check_alpha_equivalence(α1, α2)
    @test is_equiv == true
    @test c == 3
end

@testset "Search Integration" begin
    # 测试无需 optimizer
    results, failed = search_gadgets(
        AlphaTensorMode,
        loader,
        constraints
        # 注意：没有 optimizer 参数！
    )
    
    @test results isa Vector
end
```

---

## 🔍 关键实现细节

### 1. Bitmask 编码

```julia
# 顶点配置用 UInt32 表示
state = UInt32(0b101)  # 顶点 1 和 3 被选中

# 检查第 i 个顶点（0-based）
is_selected = (state >> i) & 0x1 == 1

# Pin 配置编码
# 假设 pins = [2, 5, 7]
# pin_config 的第 i 位对应 pins[i+1]
pin_config = UInt32(0b011)  # pins[1]=0, pins[2]=1, pins[3]=1
# 即：顶点2不选，顶点5选，顶点7选
```

### 2. MIS 枚举

```julia
# GadgetSearch.jl 已有实现
mis_masks, count = find_maximal_independent_sets(graph)
# 返回：
# - mis_masks: Vector{UInt32}，每个元素是一个 MIS 的 bitmask
# - count: MIS 的数量

# 计算 MIS 的大小（基数）
for mask in mis_masks
    size = count_ones(mask)
end
```

### 3. 子图诱导

```julia
# 创建仅包含特定顶点的子图
vertices = [1, 3, 5, 7]
subgraph, vmap = induced_subgraph(graph, vertices)

# vmap[i] 是子图顶点 i 在原图中的编号
```

### 4. 错误处理

```julia
# 边界情况
if isempty(feasible_interior)
    return 0  # 没有可行的内部顶点
end

if isempty(target_states)
    return nothing  # 无效输入
end

# Pin 数量限制
if length(pins) > 6
    @warn "Large number of pins ($(length(pins))), may be slow"
end
```

---

## ✅ 验证清单

实现完成后，检查：

- [ ] `compute_reduced_alpha_tensor` 对简单图（路径、三角形）返回正确值
- [ ] `check_alpha_equivalence` 正确识别等价和非等价的 α-tensor
- [ ] `infer_pattern_alpha` 从目标状态正确推断 pattern
- [ ] `verify_gadget_via_alpha_tensor` 返回 uniform weights
- [ ] `AlphaTensorMode` 在 `search_gadgets` 中工作
- [ ] 不需要提供 `optimizer` 参数
- [ ] 所有测试通过
- [ ] 与 `RydbergUnweightedModel` 找到相似的 gadgets（对简单情况）

---

## 📚 参考实现位置

如果需要查看完整实现：

```
GadgetSearch.jl/
├── src/
│   ├── core/
│   │   ├── alpha_tensor.jl      # 核心算法（~350 行）
│   │   └── search.jl            # 集成（修改 ~100 行）
│   └── GadgetSearch.jl          # 导出（修改 ~10 行）
├── test/
│   └── core/
│       └── alpha_tensor_test.jl # 测试（~250 行）
└── examples/
    └── alpha_tensor_example.jl  # 示例（~200 行）
```

---

## 🎯 最终 API

### 用户视角

```julia
using GadgetSearch

# 方式 1: 搜索 gadgets（自动推断 pattern）
results, failed = search_gadgets(
    AlphaTensorMode,
    loader,
    [TruthTableConstraint(BitMatrix([...]))]
)

# 方式 2: 直接计算 α-tensor
α = compute_reduced_alpha_tensor(graph, pins)

# 方式 3: 检查两个 gadget 的等价性
is_equiv, overhead = check_alpha_equivalence(α1, α2)
```

### 与其他模式对比

```julia
# 需要 optimizer
using HiGHS
search_gadgets(RydbergModel, loader, constraints; optimizer=HiGHS.Optimizer)

# 不需要 optimizer（简单检查）
search_gadgets(RydbergUnweightedModel, loader, constraints)

# 不需要 optimizer（α-tensor 检查，更严格）
search_gadgets(AlphaTensorMode, loader, constraints)
```

---

## ⚠️ 注意事项

1. **Pin 数量限制**：建议 ≤ 4，最多 6（因为复杂度 O(2^|pins|)）
2. **图的大小**：最多 32 个顶点（UInt32 限制）
3. **数据集**：必须是 King's lattice/triangular lattice 的子图
4. **等价性**：只检查常数差异，不检查其他性质

---

## 📊 预期性能

对于典型的 gadget 搜索（3 pins，10 顶点图）：

- `compute_reduced_alpha_tensor`: ~1ms
- `check_alpha_equivalence`: <0.1ms
- 整体搜索速度：与 `RydbergUnweightedModel` 相当
- 比 `RydbergModel`（需要优化器）快 10-100 倍

---

## 🔗 总结

| 问题 | 答案 | 实现方式 |
|------|------|----------|
| TruthTable ↔ α-tensor | 互补概念 | `infer_pattern_alpha` 转换 |
| α-tensor 计算 | 暴力枚举 | 直接用 `find_maximal_independent_sets` |
| Pattern 定义 | 动态推断 | 从 target_states 推断 |
| King's subgraph | 限制数据集 | 用 `generate_full_grid_udg` |
| 架构集成 | Dispatch | `AlphaTensorMode <: EnergyModel` |

**关键原则**：复用现有代码，保持 API 一致，用 dispatch 实现多态。

