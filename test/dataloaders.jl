@testset "dataloaders" begin
    file_path = pkgdir(GadgetSearch, "datasets", "logic_gates", "2in2out.json")
    res_dict = load_gadget(file_path)
    check_gadget([2,2], res_dict[0])

    gadget_path = pkgdir(GadgetSearch, "datasets", "udgs", "2in2out.json")
    res_dict = load_grid_gadget(gadget_path)

    save_path = pkgdir(GadgetSearch, "test", "test.pdf")
    plot_single_gadget(res_dict[1], save_path)
    @test isfile(save_path)
    rm(save_path)
end