@testset "assign.jl" begin
    
    @testset "assign_secondary_structure!" begin
        protein = pdb_to_protein("data/1ASS.pdb")
        @test !has_assigned_ss(protein)
        assign_secondary_structure!(protein)
        @test has_assigned_ss(protein)
    end

    @testset "assign_secondary_structure" begin
        protein = pdb_to_protein("data/1ASS.pdb")
        @test !has_assigned_ss(protein)
        new_protein = assign_secondary_structure(protein)
        @test !has_assigned_ss(protein)
        @test has_assigned_ss(new_protein)
    end

end