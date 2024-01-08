export Frames

using LinearAlgebra
using Statistics: mean
using Rotations: QuatRotation, params

"""
    Frames{T <: Real} <: AbstractVector{Tuple{QuatRotation{T}, Vector{T}}}

The `Frames` type is designed to efficiently store and manipulate the rotation and translation of a set of frames.
"""
struct Frames{T <: Real} <: AbstractVector{Tuple{QuatRotation{T}, Vector{T}}}
    rotations::Matrix{T}
    locations::Matrix{T}

    function Frames{T}(rotations::Matrix{T}, locations::Matrix{T}) where T <: Real
        size(rotations, 2) == size(locations, 2) || throw(ArgumentError("rotations and locations must have the same number of columns"))
        size(rotations, 1) == 4 || throw(ArgumentError("rotations must be a 4xN quaternion matrix"))
        size(locations, 1) == 3 || throw(ArgumentError("locations must be a 3xN 3D coordinates matrix"))
        new{T}(rotations, locations)
    end
end

Frames(rotations::Matrix{T}, locations::Matrix{T}) where T <: Real = Frames{T}(rotations, locations)

Base.length(frames::Frames) = size(frames.rotations, 2)
Base.size(frames::Frames) = Tuple(length(frames))
Base.getindex(frames::Frames, i::Integer) = QuatRotation(frames.rotations[:, i]), frames.locations[:, i]

centroid(P::AbstractMatrix{<:Real}) = mean(P, dims=2)

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

function Frames(backbone::Backbone{T}, ideal_coords::AbstractMatrix{<:Real}) where T <: Real
    ideal_coords = T.(ideal_coords)
    num_frame_points = size(ideal_coords, 2)
    L, r = divrem(length(backbone), num_frame_points)
    iszero(r) || throw(ArgumentError("Backbone length must be a multiple of the number of points in a frame ($num_frame_points)"))
    rotations = Matrix{T}(undef, 4, L)
    locations = Matrix{T}(undef, 3, L)
    all_raw_coords = reshape(backbone.coords, 3, num_frame_points, L)
    for (i, raw_coords) in enumerate(eachslice(all_raw_coords, dims=3))
        rotmat, _, raw_centroid = kabsch_algorithm(ideal_coords, raw_coords)
        rotations[:, i] = params(QuatRotation(rotmat))
        locations[:, i] = raw_centroid
    end
    return Frames(rotations, locations)
end

function Backbone(frames::Frames{T}, ideal_coords::AbstractMatrix{<:Real}) where T <: Real
    ideal_coords = T.(ideal_coords)
    num_frame_points = size(ideal_coords, 2)
    L = length(frames)
    approx_raw_coords = Array{T}(undef, 3, num_frame_points, L)
    ideal_centroid = centroid(ideal_coords)
    for (i, (rot, raw_centroid)) in enumerate(zip(eachcol(frames.rotations), eachcol(frames.locations)))
        quatrot = QuatRotation(rot)
        approx_raw_coords[:, :, i] = quatrot * (ideal_coords .- ideal_centroid) .+ raw_centroid
    end
    return Backbone(approx_raw_coords)
end
