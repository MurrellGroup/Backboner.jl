@testset "protein.jl" begin

    @testset "Protein" begin
        A = ProteinChain("A", Backbone(randn(3, 4, 5)))
        B = ProteinChain("B", Backbone(randn(3, 4, 6)))
        protein = [A, B]
        @test protein[1] == protein["A"] == A
        @test protein[2] == protein["B"] == B
        @test length.(protein) == [5, 6]
        @test !has_assigned_ss(protein)
    end

end