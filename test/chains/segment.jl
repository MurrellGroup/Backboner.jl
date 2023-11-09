@testset "segment.jl" begin

    @testset "Segment" begin
        chain = Chain("A", randn(3, 3, 3))
        chain.ssvector .= [Loop, Helix, Strand]
        segment = Segment{Loop}(chain, 1:1)
    end

    @testset "extend_segment" begin
        chain = Chain("A", randn(3, 3, 5))
        chain.ssvector .= [Loop, Helix, Strand, Strand, Helix]
        segment = Segment{Strand}(chain, 3:4)
        @test extend_segment(segment, 0:3) == Segment{MiSSing}(chain, 2:5)
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