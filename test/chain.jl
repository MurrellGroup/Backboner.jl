@testset "chain" begin

    @testset "chain.jl" begin

        coords = randn(3, 4, 5)
        backbone = Backbone(coords)
        chain = Chain("A", backbone)
        @test chain isa Chain{Float64}
        @test chain.id == "A"
        @test chain.backbone.coords == coords
        @test chain.ssvector == fill(Unassigned, length(chain))
        @test has_missing_ss(chain)
        @test length(chain) == 5
        @test size(chain) == (5,)
        @test Chain(remove_column(backbone, 4)).backbone == add_oxygens(remove_column(backbone, 4))
        @test Chain(backbone).id == ""

        @test summary(chain) == "Chain A with 5 residues"

        io = IOBuffer()
        show(io, chain)
        @test String(take!(io)) == summary(chain)

    end

end