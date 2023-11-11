@testset "rotations.jl" begin

    @testset "Backbone constructors" begin
        locations = Float64[0; 0; 0;;]
        
        # 90 degree rotation around x axis
        quatrots = [cos(π/2/2); sin(π/2/2); 0; 0;;]
        rot_matrices = Float64[1 0 0; 0 0 -1; 0 1 0;;;]

        @test all(isapprox.(Backbone(locations, quatrots), Backbone(locations, rot_matrices), atol=1e-6))
    end

    @testset "locs_and_rots" begin
        protein = pdb_to_protein("data/1ASS.pdb")
        backbone = remove_column(protein[1].backbone, 4)
        @test all(Backbone(locs_and_rots(backbone)...) .≈ backbone)
    end

end