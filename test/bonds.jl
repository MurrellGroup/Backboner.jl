@testset "bonds.jl" begin

    @testset "vectors and lengths" begin

        backbone = Backbone([
            0.0 1.0 2.0 5.0
            0.0 2.0 4.0 8.0
            0.0 2.0 4.0 4.0
        ])

        @testset "get_atom_displacements" begin
            @test Backboner.get_atom_displacements(backbone, 1, 1, 2) == [1.0 3.0; 2.0 4.0; 2.0 0.0]
            @test Backboner.get_atom_displacements(backbone, 2, 1, 2) == [1.0; 2.0; 2.0;;]
        end

        @testset "get_atom_distances" begin
            @test Backboner.get_atom_distances(backbone, 1, 1, 2) == [3.0 5.0]
            @test Backboner.get_atom_distances(backbone, 2, 1, 2) == [3.0;;]
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

        @testset "dihedrals" begin
            @test get_dihedrals(backbone) ≈ [π/2, π, -π/2, π/4]
        end
    
    end

    @testset "ChainedBonds" begin
        protein = Backboner.Protein.readpdb("data/1ASS.pdb")
        chain = protein["A"]
        backbone = chain.backbone
        bonds = ChainedBonds(backbone)

        @testset "invertibility" begin
            bonds1 = ChainedBonds(backbone)
            bonds2 = ChainedBonds(Backbone(ChainedBonds(backbone)))
            @test bonds1.lengths ≈ bonds2.lengths && bonds1.angles ≈ bonds2.angles && bonds1.dihedrals ≈ bonds2.dihedrals
        end

        io = IOBuffer()
        show(io, bonds)
        @test String(take!(io)) == "ChainedBonds{Float64, Vector{Float64}} with 455 bonds, 454 angles, and 453 dihedrals"
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
        new_backbone = append_bonds(subbackbone, bonds.lengths[end-k+1:end], bonds.angles[end-k+1:end], bonds.dihedrals[end-k+1:end])
        @test new_backbone.coords ≈ backbone.coords
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
        new_backbone = prepend_bonds(subbackbone, bonds.lengths[1:n-k], bonds.angles[1:n-k], bonds.dihedrals[1:n-k])
        @test new_backbone.coords ≈ backbone.coords
    end

end