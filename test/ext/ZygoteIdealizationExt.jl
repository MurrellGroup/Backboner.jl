import Zygote

@testset "ZygoteIdealizationExt" begin

    @testset "idealize" begin
        backbone = Protein.readpdb("data/1ZAK.pdb")["A"].backbone
        idealized_backbone = idealize(backbone, Float64[1.46, 1.52, 1.33], Float64[1.94, 2.04, 2.13])
        idealized_bonds = ChainedBonds(idealized_backbone)
        @test all(isapprox.(get_bond_lengths(idealized_bonds), [ideal for (_, ideal) in zip(1:length(backbone)-1, Iterators.cycle(Float64[1.46, 1.52, 1.33]))], atol=1e-3))
        @test all(isapprox.(get_bond_angles(idealized_bonds), [ideal for (_, ideal) in zip(1:length(backbone)-2, Iterators.cycle(Float64[1.94, 2.04, 2.13]))], atol=1e-3))
        @test backbone == idealize(backbone, Float64[0], Float64[0], mask_tolerance=0.5)
    end

end