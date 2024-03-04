import Rotations: QuatRotation, params

@testset "frames.jl" begin

    frames = Frames(
        [
            params(QuatRotation([0 0 1; 1 0 0; 0 1 0])^1);;
            params(QuatRotation([0 0 1; 1 0 0; 0 1 0])^2);;
            params(QuatRotation([0 0 1; 1 0 0; 0 1 0])^3);;
        ],
        [
            0 10 100;
            0 0 0;
            0 0 0;
        ]
    )

    @test frames[1] == Frame{Float64}(params(QuatRotation([0 0 1; 1 0 0; 0 1 0])), [0.0, 0.0, 0.0])
    @test frames[1] == Frame{Float64}(params(QuatRotation([0 0 1; 1 0 0; 0 1 0])), [0, 0, 0])
    @test frames[1] == Frame(params(QuatRotation([0 0 1; 1 0 0; 0 1 0])), [0, 0, 0])

    @test length(frames) == 3
    @test size(frames) == (3,)
    @test frames[1] == Frame{Float64}(QuatRotation([0 0 1; 1 0 0; 0 1 0]), [0.0, 0.0, 0.0])

    standard_coords = [3 1 -4; 1 -1 0; 0 0 0]
    backbone = Backbone(frames, standard_coords)

    @test backbone.coords == [
        0.0   0.0   0.0  11.0  9.0  10.0  103.0  101.0  96.0
        3.0   1.0  -4.0   0.0  0.0   0.0    1.0   -1.0   0.0
        1.0  -1.0   0.0   3.0  1.0  -4.0    0.0    0.0   0.0
    ]

    @test Frames(backbone, standard_coords) != frames # due to numerical error
    @test all(isapprox(f1.rotation, f2.rotation; atol=1e-10) && isapprox(f1.location, f2.location; atol=1e-10) for (f1, f2) in zip(Frames(backbone, standard_coords), frames))
    @test Frames(frames.rotations, frames.locations) == frames

    @testset "constructor with rotmats" begin
        rotations = frames.rotations
        rotmats = stack(collect(QuatRotation(rot)) for rot in eachcol(rotations))
        locations = frames.locations
        @test Frames(rotmats, locations) == frames
    end

end