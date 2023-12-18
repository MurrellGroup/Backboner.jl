@testset "chain" begin

    @testset "chain.jl" begin

        coords = randn(3, 3, 5)
        backbone = Backbone(coords)
        chain = ProteinChain("A", backbone)
        @test chain.id == "A"
        @test chain.backbone.coords == coords
        @test chain.aavector == fill('G', length(chain))
        @test chain.ssvector == fill(' ', length(chain))
        @test !has_assigned_ss(chain)
        @test length(chain) == 5
        @test size(chain) == (5,)
        @test ProteinChain(backbone).id == "_"

        protein = [chain]
        @test protein["A"] == protein[:A] == protein[1]
        
        @test chain[1] == Residue(1, backbone, 'G', ' ')

        @test summary(chain) == "ProteinChain A with 5 residues"

        io = IOBuffer()
        show(io, chain)
        @test String(take!(io)) == summary(chain)

    end

    @testset "atom_distances" begin

        protein = pdb_to_protein("data/1ASS.pdb")
        backbone = protein["A"].backbone
        cn_distances = carbonyl_nitrogen_distances(backbone)
        @test length(cn_distances) == length(backbone) - 1
        @test 1.32 <= sum(cn_distances) / length(cn_distances) <= 1.34

    end

end