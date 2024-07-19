export Frames

using LinearAlgebra
using NNlib

centroid(A::AbstractArray{<:Real}; dims=2) = sum(A; dims) ./ size(A, 2)

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

Frames(rotations::A,translations::B) where {T<:Real,A<:AbstractArray{T,3},B<:AbstractMatrix{T}} = Frames{T,A,B}(rotations, translations)

function Frames(rotations::AbstractArray{<:Real,3}, translations::AbstractArray{<:Real})
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
