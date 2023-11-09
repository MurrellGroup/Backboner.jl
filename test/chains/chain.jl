@testset "chain.jl" begin

    coords = randn(3, 3, 3)
    chain = Chain("A", coords)
    @test chain isa Chain{3, Float64}
    @test chain.id == "A"
    @test chain.coords == coords
    @test chain.ssvector == fill(MiSSing, 3)
    @test chain[1:1] == Chain("A", coords[:, :, 1:1])
    @test chain[1:1].id == "A[1:1]"
    @test has_missing_ss(chain)
    @test length(chain) == 3
    @test remove_column(chain, 3).coords == coords[:, 1:2, :]

end