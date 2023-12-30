@testset "chain" begin

    @testset "chain.jl" begin

        coords = randn(3, 15)
        backbone = Backbone(coords)
        chain = Chain("A", backbone)
        @test chain.id == "A"
        @test chain.backbone.coords == coords
        @test chain.aavector == fill('G', length(chain))
        @test chain.ssvector == fill(' ', length(chain))
        @test !has_assigned_ss(chain)
        @test length(chain) == 5
        @test size(chain) == (5,)
        @test Chain(backbone).id == "_"

        protein = [chain]
        @test protein["A"] == protein[:A] == protein[1]
        
        @test chain[1] == Residue(1, 'G', ' ')

        @test summary(chain) == "Chain A with 5 residues"

        io = IOBuffer()
        show(io, chain)
        @test String(take!(io)) == summary(chain)

    end

    @testset "atom coords" begin

        protein = readpdb("data/1ASS.pdb")
        chain = protein["A"]
        L = length(chain)
        @test size(nitrogen_coords(chain), 2) == L
        @test size(alphacarbon_coords(chain), 2) == L
        @test size(carbonyl_coords(chain), 2) == L

    end

    @testset "dihedrals" begin

        protein = readpdb("data/1ASS.pdb")
        chain = protein["A"]
        L = length(chain)
        @test length(psi_angles(chain)) == L - 1
        @test length(phi_angles(chain)) == L - 1
        @test length(omega_angles(chain)) == L - 1

    end

    @testset "atom_distances" begin

        protein = readpdb("data/1ASS.pdb")
        chain = protein["A"]

        @testset "nitrogen_alphacarbon_distances" begin
            NCa_distances = nitrogen_alphacarbon_distances(chain)
            @test length(NCa_distances) == length(chain)
            @test 1.44 <= sum(NCa_distances) / length(NCa_distances) <= 1.48
        end

        @testset "alphacarbon_carbonyl_distances" begin
            CaC_distances = alphacarbon_carbonyl_distances(chain)
            @test length(CaC_distances) == length(chain)
            @test 1.50 <= sum(CaC_distances) / length(CaC_distances) <= 1.54
        end

        @testset "carbonyl_nitrogen_distances" begin
            CN_distances = carbonyl_nitrogen_distances(chain)
            @test length(CN_distances) == length(chain) - 1
            @test 1.31 <= sum(CN_distances) / length(CN_distances) <= 1.35
        end

    end

end