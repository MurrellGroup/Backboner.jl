@testset "assign.jl" begin
    
    @testset "assign_secondary_structure!" begin
        backbone = load_pdb_backbone("data/1ASS.pdb")
        assign_secondary_structure!(backbone)
        @test !any(==(MiSSing), [chain.ssvector for chain in backbone])
    end

end