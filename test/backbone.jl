@testset "backbone.jl" begin

    coords = zeros(Float64, 3, 15)
    backbone = Backbone{Float64, Matrix{Float64}}(coords)

    @test Backbone(coords) == backbone

    @test length(backbone) == 15

    @test backbone[1] == coords[:, 1]
    @test backbone[1:1] == Backbone(coords[:, 1:1])
    @test view(backbone, 2:3)[1] == backbone[2]
    backbone[2] = ones(3)
    @test backbone[2] == ones(3)

    @test hash(backbone) != hash(Backbone{Float64, Matrix{Float64}}(coords .+ 1))

end