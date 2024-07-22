using Backboner: quaternions_to_rotation_matrices, rotation_matrices_to_quaternions
using LinearAlgebra

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
        
        # Helper functions
        function quaternion_multiply(q1, q2)
            w1, x1, y1, z1 = q1
            w2, x2, y2, z2 = q2
            return [
                w1*w2 - x1*x2 - y1*y2 - z1*z2,
                w1*x2 + x1*w2 + y1*z2 - z1*y2,
                w1*y2 - x1*z2 + y1*w2 + z1*x2,
                w1*z2 + x1*y2 - y1*x2 + z1*w2
            ]
        end

        quaternion_conjugate(q) = [q[1], -q[2], -q[3], -q[4]]

        N = 10
        _Q = randn(4, N)
        Q = _Q ./ sqrt.(sum(abs2, _Q, dims=1))
        R = quaternions_to_rotation_matrices(Q)
        Q2 = rotation_matrices_to_quaternions(R)
        R2 = quaternions_to_rotation_matrices(Q2)
    
        @testset "Conversion Consistency" begin
            @test R ≈ R2
            @test all((Q .≈ Q2) .| (Q .≈ -Q2))
        end
    
        @testset "Rotation Matrix Properties" begin
            for i in 1:N
                # Test orthogonality
                @test R[:,:,i] * R[:,:,i]' ≈ I
                # Test proper rotation (determinant should be 1)
                @test det(R[:,:,i]) ≈ 1
                # Test that it actually rotates a vector
                v = randn(3)
                @test norm(R[:,:,i] * v) ≈ norm(v)
                @test R[:,:,i] * v ≉ v  # not approximately equal
            end
        end
    
        @testset "Quaternion Properties" begin
            for i in 1:N
                # Test unit norm
                @test norm(Q[:,i]) ≈ 1
                # Test rotation property
                v = randn(3)
                qv = [0; v]
                rotated_v = quaternion_multiply(quaternion_multiply(Q[:,i], qv), quaternion_conjugate(Q[:,i]))[2:4]
                @test norm(rotated_v) ≈ norm(v)
                @test rotated_v ≉ v  # not approximately equal
            end
        end
    end

end