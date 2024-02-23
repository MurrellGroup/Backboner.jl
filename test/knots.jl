@testset "knots.jl" begin
    
    coords = Float64[                                                                                                                                                                                                                    
        0 4 2  1  1  4
        0 0 2 -1 -1  2
        0 0 0 -1  1 -1
    ]

    backbone = Backbone(coords)

    @test is_knotted(backbone)
    @test !is_knotted(backbone[1:end-1])

    @test is_knotted(Protein.readpdb("data/1YVE.pdb")["I"].backbone)
    @test !is_knotted(Protein.readpdb("data/1ASS.pdb")["A"].backbone)
    @test !is_knotted(Protein.readpdb("data/1ZAK.pdb")["A"].backbone)

end