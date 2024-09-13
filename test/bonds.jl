@testset "bonds.jl" begin

    @testset "vectors and lengths" begin

        backbone = Backbone([
            0.0 1.0 2.0 5.0
            0.0 2.0 4.0 8.0
            0.0 2.0 4.0 4.0
        ])

        @testset "get_displacements" begin
            @test Backboner.get_displacements(backbone, 1, 1, 2) == [1.0 3.0; 2.0 4.0; 2.0 0.0]
            @test Backboner.get_displacements(backbone, 2, 1, 2) == [1.0; 2.0; 2.0;;]
        end

        @testset "get_distances" begin
            @test Backboner.get_distances(backbone, 1, 1, 2) == [3.0 5.0]
            @test Backboner.get_distances(backbone, 2, 1, 2) == [3.0;;]
        end

        @testset "get_bond_vectors" begin
            @test get_bond_vectors(backbone) == [
                1.0 1.0 3.0
                2.0 2.0 4.0
                2.0 2.0 0.0
            ]
        end

        @testset "get_bond_lengths" begin
            @test get_bond_lengths(backbone) == [3.0, 3.0, 5.0]
        end

    end

    @testset "angles" begin

        coords = [                     # naturally improbable edge case with two consecutive parallel vectors. might result in a NaN depending on the implementation
            0.0 1.0 1.0 1.0 1.0 2.0 2.0#  2.0;
            0.0 0.0 1.0 1.0 2.0 2.0 1.0#  0.0;
            0.0 0.0 0.0 1.0 1.0 1.0 0.0# -1.0
        ]

        backbone = Backbone(coords)

        @testset "bond angles" begin
            @test get_bond_angles(backbone) ≈ [π/2, π/2, π/2, π/2, π/2]
        end

        @testset "torsion angles" begin
            @test get_torsion_angles(backbone) ≈ [π/2, π, -π/2, π/4]
        end
    
    end

    @testset "ChainedBonds" begin
        backbone = Backbone(rand(3,100))
        bonds = ChainedBonds(backbone)

        @testset "invertibility" begin
            bonds1 = ChainedBonds(backbone)
            bonds2 = ChainedBonds(Backbone(ChainedBonds(backbone)))
            @test get_bond_lengths(bonds1) ≈ get_bond_lengths(bonds2)
            @test get_bond_angles(bonds1) ≈ get_bond_angles(bonds2)
            @test get_torsion_angles(bonds1) ≈ get_torsion_angles(bonds2)
        end

        @test sprint(show, bonds) == "ChainedBonds{Float64, Vector{Float64}} with 99 bond lengths, 98 bond angles, and 97 torsion angles"
    end

    @testset "append_bonds" begin
        backbone = Backbone([
            0.0 1.0 1.0 1.0 1.0 2.0 2.0
            0.0 0.0 1.0 1.0 2.0 2.0 1.0
            0.0 0.0 0.0 1.0 1.0 1.0 0.0
        ])

        n = length(backbone)
        k = 3

        subbackbone = backbone[1:n-k]
        bonds = ChainedBonds(backbone)
        new_backbone = append_bonds(subbackbone, get_bond_lengths(bonds)[end-k+1:end], get_bond_angles(bonds)[end-k+1:end], get_torsion_angles(bonds)[end-k+1:end])
        @test coords(new_backbone) ≈ coords(backbone)
    end

    @testset "prepend_bonds" begin
        backbone = Backbone([
            0.0 1.0 1.0 1.0 1.0 2.0 2.0
            0.0 0.0 1.0 1.0 2.0 2.0 1.0
            0.0 0.0 0.0 1.0 1.0 1.0 0.0
        ])

        n = length(backbone)
        k = 3

        subbackbone = backbone[end-k+1:end]
        bonds = ChainedBonds(backbone)
        new_backbone = prepend_bonds(subbackbone, get_bond_lengths(bonds)[1:n-k], get_bond_angles(bonds)[1:n-k], get_torsion_angles(bonds)[1:n-k])
        @test new_backbone.coords ≈ backbone.coords
    end

end