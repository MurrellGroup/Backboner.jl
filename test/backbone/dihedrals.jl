@testset "dihedrals.jl" begin
    
    protein = pdb_to_protein("data/1ASS.pdb")
    chain_A = protein["A"]
    backbone = chain_A.backbone
    backbone3 = Backbone(atom_coord_matrix(backbone, 1:3))

    @testset "idealized lengths and angles" begin

        ideal_backbone_coords = idealize_lengths_angles(backbone3.coords)
        ks, ls, dihs = get_ks_ls_dihs(backbone3.coords)
        fixed_ideal = Backboner.fix_sequence_lengths_angles(ideal_backbone_coords, ks, ls)
        @test all(mapslices(norm, fixed_ideal .- backbone3.coords, dims=1) .< 0.1)

    end

end