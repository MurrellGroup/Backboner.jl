@test "coordmatrix.jl" begin

    backbone = load_pdb_backbone("data/1ASS.pdb")
    chain = backbone[1]
    @test nitrogen_coord_matrix(chain) == chain[:, 1, :]
    @test alphacarbon_coord_matrix(chain) == chain[:, 2, :]
    @test carbon_coord_matrix(chain) == chain[:, 3, :]
    @test oxygen_coord_matrix(chain) == chain[:, 4, :]

end