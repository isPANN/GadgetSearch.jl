using Test
using GadgetSearch

@testset "generic_rule - basic structure" begin
    # (1,1): 1 input bit, 1 output bit → 2 rows, 2 columns
    tt = GadgetSearch.generic_rule(0, (1, 1))
    @test tt isa BitMatrix
    @test size(tt) == (2, 2)
end

@testset "generic_rule - constant-0 gate (rule 0, 1-in 1-out)" begin
    # rule 0: both inputs map to output 0
    # row 1: input=0, output=0 → [0, 0]
    # row 2: input=1, output=0 → [1, 0]
    tt = GadgetSearch.generic_rule(0, (1, 1))
    @test tt[1, 1] == false && tt[1, 2] == false   # input 0 → output 0
    @test tt[2, 1] == true  && tt[2, 2] == false   # input 1 → output 0
end

@testset "generic_rule - NOT gate (rule 1, 1-in 1-out)" begin
    # rule 1: input 0 → output 1, input 1 → output 0
    tt = GadgetSearch.generic_rule(1, (1, 1))
    @test tt[1, 1] == false && tt[1, 2] == true    # input 0 → output 1
    @test tt[2, 1] == true  && tt[2, 2] == false   # input 1 → output 0
end

@testset "generic_rule - identity gate (rule 2, 1-in 1-out)" begin
    # rule 2: input 0 → output 0, input 1 → output 1
    tt = GadgetSearch.generic_rule(2, (1, 1))
    @test tt[1, 1] == false && tt[1, 2] == false   # input 0 → output 0
    @test tt[2, 1] == true  && tt[2, 2] == true    # input 1 → output 1
end

@testset "generic_rule - constant-1 gate (rule 3, 1-in 1-out)" begin
    # rule 3: both inputs map to output 1
    tt = GadgetSearch.generic_rule(3, (1, 1))
    @test tt[1, 1] == false && tt[1, 2] == true    # input 0 → output 1
    @test tt[2, 1] == true  && tt[2, 2] == true    # input 1 → output 1
end

@testset "generic_rule - 2-input 1-output AND gate" begin
    # (2,1): 2 input bits, 1 output bit → 4 rows, 3 columns
    # AND gate is rule 8 (binary: 1000 → output[3]=1 only for input=3)
    tt = GadgetSearch.generic_rule(8, (2, 1))
    @test size(tt) == (4, 3)
    # Row ordering: input 0,1,2,3 in order
    # input=0 (00) → output 0
    @test tt[1, 3] == false
    # input=1 (01) → output 0
    @test tt[2, 3] == false
    # input=2 (10) → output 0
    @test tt[3, 3] == false
    # input=3 (11) → output 1
    @test tt[4, 3] == true
end

@testset "generic_rule - error on invalid rule_id" begin
    # max rule for (1,1) is 4, so rule_id=4 should error
    @test_throws ErrorException GadgetSearch.generic_rule(4, (1, 1))
    @test_throws ErrorException GadgetSearch.generic_rule(-1, (1, 1))
end

@testset "generic_rule - show_info flag" begin
    # show_info=false should not throw
    @test_nowarn GadgetSearch.generic_rule(0, (1, 1); show_info=false)
    # show_info=true should not throw and should print output
    tt = GadgetSearch.generic_rule(0, (1, 1); show_info=true)
    @test tt isa BitMatrix
    @test size(tt) == (2, 2)
end

@testset "show_rule_info - basic output" begin
    tmp = tempname()
    open(tmp, "w") do io
        redirect_stdout(io) do
            GadgetSearch.show_rule_info(0, (1, 1))
        end
    end
    output = read(tmp, String)
    @test occursin("Input (Binary)", output)
    @test occursin("Output (Binary)", output)
    @test occursin("0", output)
    rm(tmp; force=true)
end

@testset "show_rule_info - 2-input 1-output" begin
    tmp = tempname()
    open(tmp, "w") do io
        redirect_stdout(io) do
            GadgetSearch.show_rule_info(8, (2, 1))
        end
    end
    output = read(tmp, String)
    @test occursin("Input (Binary)", output)
    lines = split(strip(output), '\n')
    # header + separator + 4 data lines
    @test length(lines) >= 6
    rm(tmp; force=true)
end

@testset "generic_rule - 2-input 2-output" begin
    # (2,2): 4 inputs, 4 outputs, max_gateid = 4^4 = 256
    tt = GadgetSearch.generic_rule(0, (2, 2))
    @test size(tt) == (4, 4)
    # rule 0: all outputs are 0
    @test all(tt[:, 3:4] .== false)
end
