@testset "io.jl" begin

    backbone = load_pdb("data/1ASS.pdb")
    @test length.(backbone) == [152]

end