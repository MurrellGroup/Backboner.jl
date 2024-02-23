@testset "backbone.jl" begin

    coords = zeros(3, 20)
    backbone = Backbone(coords)
    @test length(backbone) == 20
    @test size(backbone) == (3, 20)
    @test backbone[1] == coords[:, 1]
    @test backbone[1:2] == coords[:, 1:2]

    backbone[1] = ones(3)
    @test backbone[1] == ones(3)

    @test length(Backbone{Float32}(undef, 20)) == 20
    @test length(Backbone(randn(3, 3, 5))) == 15
    @test length(Backbone{Float32}()) == 0

end