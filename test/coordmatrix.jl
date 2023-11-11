@testset "coordmatrix.jl" begin

    protein = pdb_to_protein("data/1ASS.pdb")
    backbone = protein[1].backbone
    @test nitrogen_coord_matrix(backbone) == backbone[:, 1, :]
    @test alphacarbon_coord_matrix(backbone) == backbone[:, 2, :]
    @test carbon_coord_matrix(backbone) == backbone[:, 3, :]
    @test oxygen_coord_matrix(backbone) == backbone[:, 4, :]

end