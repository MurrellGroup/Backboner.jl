import Zygote

@testset "ZygoteIdealizationExt" begin

    @testset "idealize" begin
        backbone = Protein.readpdb("data/1ZAK.pdb")["A"].backbone

        idealized = idealize(backbone, Float32[1.46, 1.52, 1.33], Float32[1.94, 2.04, 2.13])
        @test all(isapprox.(ChainedBonds(idealized).lengths, [ideal for (_, ideal) in zip(1:length(backbone)-1, Iterators.cycle(Float32[1.46, 1.52, 1.33]))], atol=1e-3))
        @test all(isapprox.(ChainedBonds(idealized).angles, [ideal for (_, ideal) in zip(1:length(backbone)-2, Iterators.cycle(Float32[1.94, 2.04, 2.13]))], atol=1e-3))

        @test backbone == idealize(backbone, Float32[0], Float32[0], mask_tolerance=0.5)
    end

end