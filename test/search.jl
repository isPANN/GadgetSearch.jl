
# @testset "match_rows_by_pinset" begin
#     mis_result = BitMatrix([
#         0 0 1 1;
#         1 0 1 0;
#         1 1 0 0;
#         0 1 0 1;
#         1 0 1 1;
#     ])
#     truth_table = BitMatrix([
#         1 0;
#         0 1;
#         1 1;
#     ])
#     pin_set = [1, 2]

#     matched = GadgetSearch.match_rows_by_pinset(mis_result, truth_table, pin_set)
#     @test matched == [
#         [2, 5],  # rows matching [1, 0]
#         [4],     # rows matching [0, 1]
#         [3],     # rows matching [1, 1]
#     ]
# end

@testset "match_rows_by_pinset" begin
    g6string = "I??FEb_z_"
    g = GadgetSearch._parse_g6_string(g6string, BitVector(undef, 1024))
    mis_result, _ = GadgetSearch.find_maximal_independent_sets(g)
    @test length(mis_result) == 6

    tt = BitMatrix([0 0 1 1;0 1 1 1; 1 0 0 0])
    matched = GadgetSearch.match_rows_by_pinset(mis_result, tt, [3,4,5,6])
    @test matched == [[4],[1],[]]
end