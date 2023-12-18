@testset "bonds.jl" begin

    @testset "vectors and lengths" begin

        backbone = Backbone{2}([
            0.0 1.0;
            0.0 2.0;
            0.0 2.0;;;

            2.0 5.0;
            4.0 8.0;
            4.0 4.0
        ])

        @testset "get_atom_displacements" begin
            @test get_atom_displacements(backbone, 1, 2, 0) == [1.0 3.0; 2.0 4.0; 2.0 0.0]
            @test get_atom_displacements(backbone, 2, 1, 1) == [1.0; 2.0; 2.0;;]
        end

        @testset "get_atom_distance" begin
            @test get_atom_distances(backbone, 1, 2, 0) == [3.0, 5.0]
            @test get_atom_distances(backbone, 2, 1, 1) == [3.0]
        end

        @testset "get_bond_vectors" begin
            @test get_bond_vectors(backbone) == [
                1.0 1.0 3.0;
                2.0 2.0 4.0;
                2.0 2.0 0.0
            ]
        end

        @testset "get_bond_lengths" begin
            @test get_bond_lengths(backbone) == [3.0, 3.0, 5.0]
        end

    end

    @testset "angles" begin

        coords = [
            0.0 1.0 1.0 1.0 1.0 2.0 2.0  2.0;
            0.0 0.0 1.0 1.0 2.0 2.0 1.0  0.0;
            0.0 0.0 0.0 1.0 1.0 1.0 0.0 -1.0
        ]

        backbone = Backbone(reshape(coords, 3, 1, :))

        @testset "bond angles" begin
            @test get_bond_angles(backbone) ≈ [π/2, π/2, π/2, π/2, π/2, π]
        end

        @testset "dihedrals" begin
            @test get_dihedrals(backbone) ≈ [π/2, π, -π/2, π/4, 0]
        end
    
    end

    @testset "ChainedBonds" begin

        @testset "invertibility" begin
            backbone = pdb_to_protein("data/1ASS.pdb")["A"].backbone
            @test ChainedBonds(Backbone{3}(ChainedBonds(backbone))) ≈ ChainedBonds(backbone)
            @test ChainedBonds(Backbone(ChainedBonds(backbone))) ≈ ChainedBonds(backbone)
        end

    end

end