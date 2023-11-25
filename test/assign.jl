@testset "assign.jl" begin
    
    @testset "assign_secondary_structure!" begin
        protein = pdb_to_protein("data/1ASS.pdb")
        @test has_missing_ss(protein)
        assign_secondary_structure!(protein)
        @test !has_missing_ss(protein)
    end

    @testset "assign_secondary_structure" begin
        protein = pdb_to_protein("data/1ASS.pdb")
        @test has_missing_ss(protein)
        new_protein = assign_secondary_structure(protein)
        @test has_missing_ss(protein)
        @test !has_missing_ss(new_protein)
    end

end