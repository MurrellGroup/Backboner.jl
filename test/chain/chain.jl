@testset "chain" begin

    @testset "chain.jl" begin

        coords = randn(3, 4, 5)
        backbone = Backbone(coords)
        chain = Chain("A", backbone)
        @test chain isa Chain{Float64}
        @test chain.id == "A"
        @test chain.backbone.coords == coords
        @test chain.ssvector == fill(MiSSing, length(chain))
        @test has_missing_ss(chain)
        @test length(chain) == 5
        @test Chain(remove_column(backbone, 4)).backbone == add_oxygen_slice(remove_column(backbone, 4))
        @test Chain(backbone).id == ""

    end

    include("segment.jl")

end