@testset "chain" begin

    @testset "chain.jl" begin

        coords = randn(3, 4, 5)
        backbone = Backbone(coords)
        chain = ProteinChain("A", backbone)
        @test chain.id == "A"
        @test chain.backbone.coords == coords
        @test chain.aavector == fill('G', length(chain))
        @test chain.ssvector == fill(' ', length(chain))
        @test !has_assigned_ss(chain)
        @test length(chain) == 5
        @test size(chain) == (5,)
        @test ProteinChain(remove_column(backbone, 4)).backbone == add_oxygens(remove_column(backbone, 4))
        @test ProteinChain(backbone).id == "_"
        
        @test chain[1] == Residue(1, backbone, 'G', ' ')

        @test summary(chain) == "ProteinChain A with 5 residues"

        io = IOBuffer()
        show(io, chain)
        @test String(take!(io)) == summary(chain)

    end

end