@testset "chains" begin

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

    @testset "segment.jl" begin

        @testset "Segment" begin
            chain = Chain("A", randn(3, 3, 3))
            chain.ssvector .= [Loop, Helix, Strand]
            segment = Segment{Loop}(chain, 1:1)
        end

        @testset "segments" begin
            coords = randn(3, 3, 3)
            chain = Chain("B", coords)
            @test_throws ErrorException segments(chain)
            chain.ssvector .= [Loop, Helix, Helix]
            chain_segments = segments(chain)
            @test length(chain_segments) == 2
            @test chain_segments[1] isa Segment{Loop, 3, Float64}
            @test chain_segments[1].chain == chain
            @test chain_segments[1].range == 1:1
            @test chain_segments[1].coords == coords[:, :, 1:1]
            @test chain_segments[2].coords == coords[:, :, 2:3]
        end

    end

end