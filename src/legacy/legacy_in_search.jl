# function generate_pin_variants(
#     g::SimpleGraph{Int},
#     pos::Vector{Tuple{Float64, Float64}},
#     pin_set::Vector{Int},
#     bit_num::Union{Int, Tuple{Int, Int}}
# )
#     if isa(bit_num, Int)
#         return collect(Combinatorics.permutations(pin_set, bit_num))
#     else
#         total = Vector{Vector{Int}}()

#         for comb in Combinatorics.permutations(pin_set, bit_num[1]+bit_num[2])
#             has_conflict = any(has_edge(g, u, v) for (u, v) in combinations(comb, 2))
#             # @show comb
#             if !has_conflict && length(comb) ≥ 4
#                 # 判断线段是否相交：1-3 和 2-4
#                 p1, p2 = pos[comb[1]], pos[comb[3]]
#                 q1, q2 = pos[comb[2]], pos[comb[4]]
#                 if segments_intersect(p1, p2, q1, q2)
#                     push!(total, comb)
#                 end
#             end
#         end

#         return total
#         # input_len, output_len = bit_num

#         # for input in Combinatorics.combinations(pin_set, input_len)
#         #     rest = setdiff(pin_set, input)
#         #     if length(rest) ≥ output_len
#         #         for output in permutations(rest, output_len)
#         #             combined = vcat(input, output)

#         #             has_conflict = any(has_edge(g, u, v) for (u, v) in combinations(combined, 2))

#         #             if !has_conflict && length(combined) ≥ 4
#         #                 # # 判断线段是否相交：1-3 和 2-4
#         #                 # p1, p2 = pos[combined[1]], pos[combined[3]]
#         #                 # q1, q2 = pos[combined[2]], pos[combined[4]]
#         #                 # if segments_intersect(p1, p2, q1, q2)
#         #                 #     push!(total, combined)
#         #                 # end
#         #                 push!(total, combined)
#         #             end
#                 # end
#             # end
#         # end
#         return total
#     end
# end

# function segments_intersect(p1, p2, q1, q2)::Bool
#     # 向量叉积判断顺时针 / 逆时针
#     function ccw(a, b, c)
#         return (c[2] - a[2]) * (b[1] - a[1]) > (b[2] - a[2]) * (c[1] - a[1])
#     end
#     return (ccw(p1, q1, q2) != ccw(p2, q1, q2)) && (ccw(p1, p2, q1) != ccw(p1, p2, q2))
# end