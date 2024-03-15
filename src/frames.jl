export Frame, Frames

using LinearAlgebra
using Rotations: QuatRotation, params

centroid(P::AbstractMatrix{<:Real}) = vec(sum(P, dims=2)) ./ size(P, 2)

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

"""
    Frame{T <: Real}

A `Frame` is a combination of a rotation and a translation, which can be applied to a set of coordinates.
"""
struct Frame{T <: Real}
    rotation::QuatRotation{T}
    location::AbstractVector{T}

    function Frame{T}(rotation::QuatRotation{T}, location::AbstractVector{T}) where T <: Real
        length(location) == 3 || throw(ArgumentError("location must be a 3-vector"))
        return new{T}(rotation, location)
    end
end

Frame{T}(rotation::AbstractVecOrMat{T}, location::AbstractVector{T}) where T <: Real = Frame{T}(QuatRotation(rotation), location)
Frame{T}(rotation::AbstractVecOrMat{<:Real}, location::AbstractVector{<:Real}) where T <: Real = Frame{T}(T.(rotation), T.(location))
Frame(rotation::AbstractVecOrMat, location::AbstractVector) = Frame{promote_type(eltype(rotation), eltype(location))}(rotation, location)

Base.:(==)(frame1::Frame, frame2::Frame) = frame1.rotation == frame2.rotation && frame1.location == frame2.location

function (frame::Frame{T})(coords::AbstractMatrix{T}, coords_centroid::AbstractVector{T}=centroid(coords)) where T 
    return frame.rotation * (coords .- coords_centroid) .+ frame.location
end

"""
    Frames{T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{Frame{T}}

The `Frames` type is designed to efficiently store and manipulate the rotation and translation of a set of `Frame`s.
"""
struct Frames{T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{Frame{T}}
    rotations::M
    locations::M

    function Frames{T}(rotations::M, locations::M) where {T <: Real, M <: AbstractMatrix{T}}
        size(rotations, 1) == 4 || throw(ArgumentError("rotations must be a 4xN quaternion matrix"))
        size(locations, 1) == 3 || throw(ArgumentError("locations must be a 3xN 3D coordinates matrix"))
        size(rotations, 2) == size(locations, 2) || throw(ArgumentError("rotations and locations must have the same number of columns"))
        return new{T, M}(rotations, locations)
    end
end

function Frames{T}(rotmats::AbstractArray{T, 3}, locations::AbstractMatrix{T}) where T <: Real
    rotations = stack(params(QuatRotation(rotmat)) for rotmat in eachslice(rotmats, dims=3))
    return Frames{T}(rotations, locations)
end

function Frames{T}(rotations::AbstractArray{<:Real}, locations::AbstractMatrix{<:Real}) where T <: Real
    return Frames{T}(T.(rotations), T.(locations))
end

function Frames(rotations::AbstractArray{<:Real}, locations::AbstractMatrix{<:Real})
    T = promote_type(eltype(rotations), eltype(locations))
    return Frames{T}(rotations, locations)
end

Base.size(frames::Frames) = Tuple(size(frames.rotations, 2))
Base.getindex(frames::Frames{T}, i::Integer) where T = Frame{T}(QuatRotation(frames.rotations[:, i]), frames.locations[:, i])

Base.:(==)(frames1::Frames, frames2::Frames) = all(f1 == f2 for (f1, f2) in zip(frames1, frames2))

function (frames::Frames{T})(coords::AbstractMatrix{T}) where T <: Real
    coords = T.(coords)
    coords_centroid = centroid(coords)
    return stack((f -> f(coords, coords_centroid)).(frames))
end

function Frames(backbone::Backbone{<:Real}, ideal_coords::AbstractMatrix{<:Real})
    T = promote_type(eltype(backbone.coords), eltype(ideal_coords))
    backbone = Backbone{T}(backbone.coords)
    ideal_coords = T.(ideal_coords)
    L, r = divrem(length(backbone), size(ideal_coords, 2))
    iszero(r) || throw(ArgumentError("backbone length ($(length(backbone))) must be divisible of the number of frame points ($(size(ideal_coords, 2)))"))
    rotmats = similar(backbone.coords, 3, 3, L)
    locations = similar(backbone.coords, 3, L)
    for (i, noisy_coords) in enumerate(eachslice(reshape(backbone.coords, 3, size(ideal_coords, 2), :), dims=3))
        rotmats[:, :, i], _, locations[:, i] = kabsch_algorithm(ideal_coords, noisy_coords)
    end
    return Frames(rotmats, locations)
end

Backbone(frames::Frames, ideal_coords::AbstractMatrix{<:Real}) = Backbone(frames(ideal_coords))
