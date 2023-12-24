@testset "backbone.jl" begin

    coords = randn(3, 20)
    backbone = Backbone(coords)
    @test length(backbone) == 20
    @test size(backbone) == (20,)
    @test backbone[1] == coords[:, 1]

end