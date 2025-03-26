@testset "dataloaders" begin
    file_path = pkgdir(GadgetSearch, "datasets", "logic_gates", "2in2outn.json")
    load_gadget(file_path)

    gadget_path = pkgdir(GadgetSearch, "datasets", "udgs", "2in2out_m3n4n.json")
    pos_path = pkgdir(GadgetSearch, "data", "grid_udgs", "m3n4pad1_min3max12_direct4.json")
    load_grid_gadget(3, 4, 1, gadget_path, pos_path)
end