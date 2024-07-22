@testset "frames.jl" begin

    @testset "Frames" begin
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

    @testset "Quaternion conversion" begin
        N = 10
        _Q = randn(4, N)
        Q = _Q ./ sqrt.(sum(abs2, _Q, dims=1))
        R = quaternions_to_rotation_matrices(Q)
        Q2 = rotation_matrices_to_quaternions(R)
        R2 = quaternions_to_rotation_matrices(Q2)
        @test R ≈ R2
        @test all((Q .≈ Q2) .| (Q .≈ -Q2))
    end

end