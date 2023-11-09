@testset "protein.jl" begin

    @testset "Protein" begin
        A = Chain("A", Backbone(randn(3, 4, 5)))
        B = Chain("B", Backbone(randn(3, 4, 6)))
        protein = Protein([A, B])
        @test protein[1] == protein["A"] == A
        @test protein[2] == protein["B"] == B
        @test length(protein) == 2
        @test length.(protein) == [5, 6]
        @test has_missing_ss(protein)
    end

end