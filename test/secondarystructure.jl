@testset "secondarystructure.jl" begin

    @test Loop == SecondaryStructure(1)
    @test Helix == SecondaryStructure(2)
    @test Strand == SecondaryStructure(3)
    @test has_missing_ss([Loop, Helix, Strand, Unassigned])
    @test !has_missing_ss([Loop, Helix, Strand])

end