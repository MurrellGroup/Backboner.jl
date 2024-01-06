@testset "bonds.jl" begin

    @testset "vectors and lengths" begin

        backbone = Backbone([
            0.0 1.0 2.0 5.0;
            0.0 2.0 4.0 8.0;
            0.0 2.0 4.0 4.0;
        ])

        @testset "get_atom_displacements" begin
            @test get_atom_displacements(backbone, 1, 1, 2) == [1.0 3.0; 2.0 4.0; 2.0 0.0]
            @test get_atom_displacements(backbone, 2, 1, 2) == [1.0; 2.0; 2.0;;]
        end

        @testset "get_atom_distances" begin
            @test get_atom_distances(backbone, 1, 1, 2) == [3.0, 5.0]
            @test get_atom_distances(backbone, 2, 1, 2) == [3.0]
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

        backbone = Backbone(coords)

        @testset "bond angles" begin
            @test get_bond_angles(backbone) ≈ [π/2, π/2, π/2, π/2, π/2, π]
        end

        @testset "dihedrals" begin
            @test get_dihedrals(backbone) ≈ [π/2, π, -π/2, π/4, 0]
        end
    
    end

    @testset "ChainedBonds" begin

        protein = Backboner.Protein.readpdb("data/1ASS.pdb")
        chain = protein["A"]
        backbone = chain.backbone
        bonds = ChainedBonds(backbone)
        @test size(bonds) == (length(backbone) - 1,)

        @testset "invertibility" begin
            @test ChainedBonds(Backbone(ChainedBonds(backbone))) ≈ ChainedBonds(backbone)
        end

        @test ChainedBonds(@view(bonds.lengths[1:end]), bonds.angles, bonds.dihedrals) == bonds

        io = IOBuffer()
        show(io, bonds)
        @test String(take!(io)) == "ChainedBonds{Float32} with 455 bonds, 454 angles, and 453 dihedrals"

    end

end