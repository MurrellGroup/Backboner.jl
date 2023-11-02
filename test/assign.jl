@testset "assign.jl" begin
    
    @testset "assign_secondary_structure!" begin
        backbone = load_pdb_backbone("data/1ASS.pdb")
        @test has_missing_ss(backbone)
        assign_secondary_structure!(backbone)
        @test !has_missing_ss(backbone)
    end

end