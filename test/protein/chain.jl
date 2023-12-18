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
        chain = protein["A"]

        @testset "nitrogen_alphacarbon_distances" begin
            NCa_distances = nitrogen_alphacarbon_distances(chain)
            @test length(NCa_distances) == length(chain)
            @test 1.45 <= sum(NCa_distances) / length(NCa_distances) <= 1.47
        end

        @testset "alphacarbon_carbonyl_distances" begin
            CaC_distances = alphacarbon_carbonyl_distances(chain)
            @test length(CaC_distances) == length(chain)
            @test 1.51 <= sum(CaC_distances) / length(CaC_distances) <= 1.53
        end

        @testset "carbonyl_nitrogen_distances" begin
            CN_distances = carbonyl_nitrogen_distances(chain)
            @test length(CN_distances) == length(chain) - 1
            @test 1.32 <= sum(CN_distances) / length(CN_distances) <= 1.34
        end

    end

end