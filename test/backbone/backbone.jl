@testset "backbone" begin
    
    @testset "backbone.jl" begin

        coords = randn(3, 4, 5)
        backbone = Backbone(coords)
        @test backbone isa Backbone{4}
        @test size(backbone) == (3, 4, 5)
        @test length(backbone) == 20
        @test backbone[1] == coords[:, :, 1]

    end

    include("rotations.jl")
    include("bonds.jl")

end