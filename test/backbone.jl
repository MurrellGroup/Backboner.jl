@testset "backbone.jl" begin

    coords = zeros(Float64, 3, 15)
    backbone = Backbone{Float64, Matrix{Float64}}(coords)

    @test Backbone{Float64}(coords) == backbone
    @test Backbone{Float64}(reshape(coords, 3, 3, :)) == backbone
    @test Backbone{Float32}(coords) == backbone
    @test Backbone{Float32}(backbone) == backbone
    @test Backbone(coords) == backbone

    @test length(Backbone{Float64}(undef, 15)) == 15
    @test length(Backbone(randn(3, 3, 5))) == 15
    @test length(Backbone{Float64}()) == 0

    @test size(backbone) == (3, 15)
    @test backbone[1, 1] == 0.0
    @test backbone[:, 1:2] == coords[:, 1:2]
    @test view(backbone, :, 1:2) == coords[:, 1:2]
    backbone[:, 1] = ones(3)
    @test backbone[:, 1] == ones(3)

    @test length(backbone) == 15
    @test backbone[1] == coords[:, 1]
    @test backbone[1:2] == coords[:, 1:2]
    @test view(backbone, 1) == coords[:, 1]
    @test view(backbone, 1:2) == coords[:, 1:2]
    backbone[2] = ones(3)
    @test backbone[2] == ones(3)

    @test backbone == coords
    @test coords == backbone

    @test copy(backbone) == coords

    @test hash(backbone) != hash(Backbone{Float64, Matrix{Float64}}(coords .+ 1))

end