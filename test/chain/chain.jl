@testset "chain" begin

    @testset "chain.jl" begin

        coords = randn(3, 4, 5)
        chain = Chain("A", Backbone(coords))
        @test chain isa Chain{Float64}
        @test chain.id == "A"
        @test chain.backbone.coords == coords
        @test chain.ssvector == fill(MiSSing, length(chain))
        @test has_missing_ss(chain)
        @test length(chain) == 5

    end

    include("segment.jl")

end