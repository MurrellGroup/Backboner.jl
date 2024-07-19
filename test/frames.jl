@testset "frames.jl" begin

    frames = Frames(
        [
            [0 0 1; 1 0 0; 0 1 0];;;
            [0 1 0; 0 0 1; 1 0 0];;;
            [1 0 0; 0 1 0; 0 0 1];;;
        ],
        [
            0 10 100;
            0 0 0;
            0 0 0;
        ]
    )

    standard_coords = [3 1 -4; 1 -1 0; 0 0 0]
    backbone = Backbone(frames, standard_coords)

    @test backbone.coords == [
        0.0   0.0   0.0  11.0  9.0  10.0  103.0  101.0  96.0
        3.0   1.0  -4.0   0.0  0.0   0.0    1.0   -1.0   0.0
        1.0  -1.0   0.0   3.0  1.0  -4.0    0.0    0.0   0.0
    ]

    new_frames = Frames(backbone, standard_coords)
    @test isapprox(frames.rotations, new_frames.rotations; atol=1e-10)
    @test isapprox(frames.translations, new_frames.translations; atol=1e-10)

end