@testset "backbone.jl" begin
    
    backbone = Backbone(randn(3, 4, 5))
    @test size(backbone) == (3, 4, 5)
    @test length(backbone) == 5
    @test remove_column(backbone, 4).coords == backbone.coords[:, 1:3, :]

end