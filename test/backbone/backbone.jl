@testset "backbone" begin
    
    @testset "backbone.jl" begin

        coords = randn(3, 4, 5)
        backbone = Backbone(coords)
        @test size(backbone) == (3, 4, 5)
        @test length(backbone) == 5
        @test backbone[1] == coords[:, :, 1]

    end

    include("rotations.jl")
    include("bonds.jl")
    include("dihedrals.jl")

end