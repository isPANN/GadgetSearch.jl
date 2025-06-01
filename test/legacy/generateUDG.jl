@testset "generateUDG.jl" begin
    save_path = pkgdir(GadgetSearch, "test")
    generate_grid_udgs(3,3,1; save_path=save_path)
    g6_file_name = pkgdir(GadgetSearch, "test", "m3n3pad1_min3max9_direct4.g6")
    pos_file_name = pkgdir(GadgetSearch, "test", "m3n3pad1_min3max9_direct4.json")
    @test isfile(g6_file_name)
    @test isfile(pos_file_name)
    rm(g6_file_name)
    rm(pos_file_name)
end