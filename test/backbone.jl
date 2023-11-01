@testset "backbone.jl" begin

    @testset "Backbone" begin
        A = Chain("A", randn(3, 3, 3))
        B = Chain("B", randn(3, 3, 4))
        backbone = Backbone([A, B])
        @test backbone[1] == A
        @test backbone[2] == B
        @test length(backbone) == 2
        @test length.(backbone) == [3, 4]
    end

end