@testset "rotations.jl" begin

    protein = pdb_to_protein("data/1ASS.pdb")
    backbone = remove_column(protein[1].backbone, 4)
    @test all(Backbone(locs_and_rots(backbone)...) .≈ backbone)

end