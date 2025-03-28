@testset "dataloaders" begin
    file_path = pkgdir(GadgetSearch, "datasets", "logic_gates", "2in2out.json")
    res_dict = load_gadget(file_path)

    gadget_path = pkgdir(GadgetSearch, "datasets", "udgs", "2in2out_m3n4.json")
    pos_path = pkgdir(GadgetSearch, "data", "grid_udgs", "m3n4pad1_min3max12_direct4.json")
    res_dict = load_grid_gadget(3, 4, 1, gadget_path, pos_path)

    save_path = pkgdir(GadgetSearch, "test", "test.pdf")
    plot_single_gadget(res_dict[1], save_path)
    @test isfile(save_path)
    rm(save_path)
end