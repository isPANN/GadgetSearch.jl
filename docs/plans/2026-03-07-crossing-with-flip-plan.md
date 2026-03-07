# Implementation Plan: Crossing Gadget Search with Logical Flip

**Design Document:** 2026-03-07-crossing-with-flip-design.md
**Estimated Complexity:** Medium (3-4 hours)

## Overview

Implement multi-target search for crossing gadgets with logical flip variants. Core strategy: generate all non-equivalent flip combinations of CROSS graph, search using existing multi-target infrastructure, visualize results with Typst.

## Tasks

### Task 1: Implement flip_alpha_tensor function
**File:** `src/utils/flip_variants.jl` (new)
**Dependencies:** None
**Estimated time:** 15 min

Create function to flip specified dimensions of a reduced alpha tensor:
```julia
function flip_alpha_tensor(tensor::Array{T,N}, pins_to_flip::Vector{Int}) where {T,N}
```

- Input: N-dimensional tensor, list of dimensions to flip (1-indexed)
- Output: New tensor with specified dimensions inverted
- Implementation: Use array indexing with `1 .- indices` for flipped dimensions

**Acceptance criteria:**
- Correctly flips single dimension
- Correctly flips multiple dimensions
- Preserves -Inf values in correct positions

---

### Task 2: Implement symmetry-aware flip variant generation
**File:** `src/utils/flip_variants.jl`
**Dependencies:** Task 1
**Estimated time:** 30 min

Create function to generate non-equivalent flip patterns:
```julia
function generate_flip_variants(graph::SimpleGraph, boundary::Vector{Int})
```

- For CROSS: manually specify non-equivalent patterns (no flip, flip-1, flip-1-2, flip-1-3, flip-all)
- Return: Vector of (graph, boundary, flip_description) tuples
- Each tuple includes modified reduced alpha tensor

**Acceptance criteria:**
- Generates 5 non-equivalent flip patterns for CROSS
- Each pattern has correct flipped alpha tensor
- Includes descriptive labels

---

### Task 3: Implement extended CROSS variants
**File:** `src/utils/flip_variants.jl`
**Dependencies:** None
**Estimated time:** 20 min

Create function to generate CROSS with inserted nodes:
```julia
function generate_extended_cross()
```

- Generate variants: 1-5-3, 2-6-4, 1-5-3 + 2-6-4
- Return: Vector of (graph, boundary) tuples
- Boundary always [1,2,3,4]

**Acceptance criteria:**
- Generates 3 extended variants
- All graphs are valid SimpleGraphs
- Boundary pins correctly identified

---

### Task 4: Combine all targets
**File:** `src/utils/flip_variants.jl`
**Dependencies:** Tasks 2, 3
**Estimated time:** 10 min

Create master function:
```julia
function generate_all_crossing_targets()
```

- Combine base CROSS flips + extended CROSS flips
- Return: Vector{Tuple{SimpleGraph{Int}, Vector{Int}}} for multi-target search

**Acceptance criteria:**
- Returns ~20 total targets (5 base + 3×5 extended)
- All targets have 4 boundary pins
- No duplicate targets

---

### Task 5: Create search script
**File:** `examples/crossing_with_flip.jl` (new)
**Dependencies:** Task 4
**Estimated time:** 30 min

Main search script:
- Load triangular UDG datasets
- Generate all crossing targets
- Call `search_multi_target_gadgets`
- Save results to JSON

**Acceptance criteria:**
- Runs without errors
- Searches m3n3 dataset successfully
- Outputs results with target indices

---

### Task 6: Implement Typst visualization
**File:** `examples/crossing_with_flip.jl`
**Dependencies:** Task 5
**Estimated time:** 45 min

Generate Typst document:
- Summary table: target index, vertices, edges, offset
- Graph visualizations using existing plot functions
- Export SVGs, embed in Typst
- Compile to PDF on desktop

**Acceptance criteria:**
- Generates valid Typst markup
- Includes all found gadgets
- PDF renders correctly on desktop

---

### Task 7: Export flip_variants functions
**File:** `src/GadgetSearch.jl`
**Dependencies:** Tasks 1-4
**Estimated time:** 5 min

Add exports:
```julia
export flip_alpha_tensor
export generate_flip_variants
export generate_extended_cross
export generate_all_crossing_targets
```

**Acceptance criteria:**
- Functions accessible from GadgetSearch module

---

### Task 8: Test and validate
**File:** `examples/crossing_with_flip.jl`
**Dependencies:** All previous
**Estimated time:** 30 min

- Run search on small dataset
- Verify results are valid gadgets
- Check visualization output
- Validate against known BATOIDEA

**Acceptance criteria:**
- Script completes successfully
- Results pass `is_gadget_replacement` checks
- PDF generated on desktop

## Implementation Order

1. Task 1 → Task 2 → Task 3 → Task 4 (target generation)
2. Task 7 (exports)
3. Task 5 (search script)
4. Task 6 (visualization)
5. Task 8 (validation)

## Risk Mitigation

- **Symmetry analysis complexity**: Start with manual enumeration for CROSS, can automate later
- **Large target count**: Use existing inf-mask optimization, should handle 20 targets efficiently
- **Typst compilation**: Reuse existing compile_typst.jl infrastructure

## Success Criteria

- [ ] All 8 tasks completed
- [ ] Search runs on m3n3 dataset
- [ ] PDF with results generated on desktop
- [ ] All found gadgets validated as correct replacements
