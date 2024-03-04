import Zygote

@testset "idealization.jl" begin

    chain = Protein.readpdb("data/1ZAK.pdb")["A"]
    idealized = idealize(chain)
    @test all(isapprox.(ChainedBonds(idealized.backbone).lengths, [ideal for (_, ideal) in zip(1:length(chain.backbone)-1, Iterators.cycle(Protein.BACKBONE_BOND_LENGTHS))], atol=1e-3))
    @test all(isapprox.(ChainedBonds(idealized.backbone).angles, [ideal for (_, ideal) in zip(1:length(chain.backbone)-2, Iterators.cycle(Protein.BACKBONE_BOND_ANGLES))], atol=1e-3))

end