export Frames

using LinearAlgebra
using NNlib

centroid(A::AbstractArray{<:Real}; dims=2) = sum(A; dims) ./ size(A, 2)

# TODO: batched version? possible? batched svd?
function kabsch_algorithm(P::AbstractMatrix{T}, Q::AbstractMatrix{T}) where T <: Real
    size(P) == size(Q) || throw(ArgumentError("P and Q must have the same size"))
    P_centroid = centroid(P)
    Q_centroid = centroid(Q)
    P_centered = P .- P_centroid
    Q_centered = Q .- Q_centroid
    H = P_centered * Q_centered'
    F = svd(H)
    U, Vt = F.U, F.Vt
    d = sign(det(U * Vt))
    Vt[end, :] .*= d
    R = (U * Vt)'
    return R, P_centroid, Q_centroid
end

struct Frames{T<:Real,A<:AbstractArray{T,3},B<:AbstractArray{T,2}}
    rotations::A
    translations::B

    function Frames{T,A,B}(rotations::A, translations::B) where {T,A,B}
        size(rotations)[1:2] == (3,3) || throw(ArgumentError("rotations must be a 3x3xL array"))
        size(translations, 1) == 3 || throw(ArgumentError("translations must be a 3xN matrix"))
        size(rotations, 3) == size(translations, 2) || throw(ArgumentError("rotations and translations must have the same number of columns"))
        return new(rotations, translations)
    end
end

Frames(rotations::A, translations::B) where {T<:Real,A<:AbstractArray{T,3},B<:AbstractArray{T,2}} = Frames{T,A,B}(rotations, translations)

function Frames(rotations::AbstractArray{<:Real}, translations::AbstractArray{<:Real})
    T = promote_type(eltype(rotations), eltype(translations))
    return Frames(T.(rotations), T.(translations))
end

function Frames(backbone::Backbone{T}, ideal_coords::AbstractMatrix{<:Real}) where T <: Real
    backbone = Backbone(backbone.coords)
    ideal_coords = T.(ideal_coords)
    L, r = divrem(length(backbone), size(ideal_coords, 2))
    iszero(r) || throw(ArgumentError("backbone length ($(length(backbone))) must be divisible of the number of frame points ($(size(ideal_coords, 2)))"))
    rotations = similar(backbone.coords, 3, 3, L)
    translations = similar(backbone.coords, 3, L)
    for (i, noisy_coords) in enumerate(eachslice(reshape(backbone.coords, 3, size(ideal_coords, 2), :), dims=3))
        rotations[:, :, i], _, translations[:, i] = kabsch_algorithm(ideal_coords, noisy_coords)
    end
    return Frames(rotations, translations)
end

(frames::Frames{T})(coords::AbstractMatrix{T}) where T<:Real = frames.rotations ⊠ (coords .- centroid(coords)) .+ reshape(frames.translations, 3, 1, :)
(frames::Frames{T})(coords::AbstractMatrix{<:Real}) where T<:Real = frames(T.(coords))

Backbone(frames::Frames, ideal_coords::AbstractMatrix{<:Real}) = Backbone(frames(ideal_coords))

### Quaternion support

# takes a batch of unit quaternions in a 4xN matrix and returns a batch of rotation matrices in a 3x3xN array
function quaternions_to_rotation_matrices(q::AbstractArray{<:Real,2})
    size(q, 1) == 4 || throw(ArgumentError("Quaternion batch must have shape 4xN"))

    sx = 2q[1, :] .* q[2, :]
    sy = 2q[1, :] .* q[3, :]
    sz = 2q[1, :] .* q[4, :]
    xx = 2q[2, :] .^ 2
    xy = 2q[2, :] .* q[3, :]
    xz = 2q[2, :] .* q[4, :]
    yy = 2q[3, :] .^ 2
    yz = 2q[3, :] .* q[4, :]
    zz = 2q[4, :] .^ 2  

    r1 = 1 .- (yy + zz)
    r2 = xy - sz
    r3 = xz + sy
    r4 = xy + sz
    r5 = 1 .- (xx + zz)
    r6 = yz - sx
    r7 = xz - sy
    r8 = yz + sx
    r9 = 1 .- (xx + yy)

    return reshape(vcat(r1', r4', r7', r2', r5', r8', r3', r6', r9'), 3, 3, :)
end

Frames(rotations::AbstractArray{T,2}, translations::AbstractArray{T,2}) = Frames(quaternions_to_rotation_matrices(rotations), translations)

# takes a batch of rotation matrices in a 3x3xN array and returns a batch of unit quaternions in a 4xN matrix
function rotation_matrices_to_quaternions(R::AbstractArray{<:Real,3})
    size(q)[1:2] == (3,3) || throw(ArgumentError("Rotation matrix batch must have shape 3x3xN"))
    # 1x1xN
    r11, r12, r13 = R[1:1, 1:1, :], R[1:1, 2:2, :], R[1:1, 3:3, :]
    r21, r22, r23 = R[2:2, 1:1, :], R[2:2, 2:2, :], R[2:2, 3:3, :]
    r31, r32, r33 = R[3:3, 1:1, :], R[3:3, 2:2, :], R[3:3, 3:3, :]

    # 4x1xN
    q0 = [1 .+ r11 + r22 + r33
          r32 - r23
          r13 - r31
          r21 - r12]
    q1 = [r32 - r23
          1 .+ r11 - r22 - r33
          r12 + r21
          r13 + r31]
    q2 = [r13 - r31
          r12 + r21
          1 .- r11 + r22 - r33
          r23 + r32]
    q3 = [r21 - r12
          r13 + r31
          r23 + r32
          1 .- r11 - r22 + r33]

    # 4x4xN
    Q = hcat(q0, q1, q2, q3)

    # 1x4xN, norm of each quaternion
    norms = sqrt.(sum(abs2, Q, dims=1))

    exp_norms = exp.(norms)
    # 1x4xN, norm weights
    weights = exp_norms ./ sum(exp_norms, dims=3)

    # batched matmul, 4x1xN
    q = Q ⊠ reshape(weights, 4, 1, :)
    q_normalized = q ./ sqrt.(sum(abs2, q, dims=1))
    
    return reshape(q_normalized, 4, :)
end