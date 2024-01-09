@testset "backbone.jl" begin

    coords = randn(3, 20)
    backbone = Backbone(coords)
    @test length(backbone) == 20
    @test size(backbone) == (20,)
    @test backbone[1] == coords[:, 1]
    @test backbone[1:2] == Backbone(coords[:, 1:2])

    backbone[1] = [0.0, 0.0, 0.0]
    backbone[1:2] = [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test backbone[1] == [0.0, 0.0, 0.0]
    @test backbone[1:2] == Backbone([0.0 0.0; 0.0 0.0; 0.0 0.0])

    @test length(Backbone{Float32}(undef, 20)) == 20
    @test length(Backbone(randn(3, 3, 5))) == 15
    @test length(Backbone{Float32}()) == 0

    io = IOBuffer()
    show(io, backbone)
    @test String(take!(io)) == "20-element Backbone{Float64, Matrix{Float64}}"


end