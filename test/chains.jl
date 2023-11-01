@testset "chains" begin

    @testset "chain.jl" begin
        coords = randn(3, 3, 3)
        chain = Chain("A", coords)
        @test chain isa Chain{3, Float64}
        @test chain.id == "A"
        @test chain.coords == coords
        @test chain.ssvector == fill(MiSSing, 3)
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
            @test_throws ErrorException chain_segments(chain)
            chain.ssvector .= [Loop, Helix, Helix]
            segments = chain_segments(chain)
            @test length(segments) == 2
            @test segments[1] isa Segment{Loop, 3, Float64}
            @test segments[1].chain == chain
            @test segments[1].range == 1:1
            @test segments[1].coords == coords[:, :, 1:1]
            @test segments[2].coords == coords[:, :, 2:3]
        end

    end

end