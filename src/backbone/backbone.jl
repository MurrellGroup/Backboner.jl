export Backbone
export get_atom
export num_atoms

"""
    Backbone{T <: Real} <: AbstractMatrix{T}

The `Backbone` type is designed to efficiently store and manipulate the three-dimensional coordinates of backbone atoms.
"""
struct Backbone{T <: Real} <: AbstractVector{AbstractVector{T}}
    coords::AbstractMatrix{T}

    function Backbone(coords::AbstractMatrix{T}) where T
        @assert size(coords, 1) == 3 "coords must have 3 coordinates per atom"
        return new{T}(coords)
    end

    function Backbone(coords::AbstractArray{T, 3}) where T
        @assert size(coords, 1) == 3 "coords must have 3 coordinates per atom"
        return new{T}(reshape(coords, 3, :))
    end
end

@inline Base.:(==)(backbone1::Backbone, backbone2::Backbone) = backbone1.coords == backbone2.coords
@inline Base.:(≈)(backbone1::Backbone, backbone2::Backbone) = backbone1.coords ≈ backbone2.coords
@inline Base.length(backbone::Backbone) = size(backbone.coords, 2)
@inline Base.size(backbone::Backbone) = Tuple(length(backbone))
@inline Base.getindex(backbone::Backbone, i::Integer) = view(backbone.coords, :, i)
@inline Base.getindex(backbone::Backbone, r::AbstractVector{<:Integer}) = Backbone(view(backbone.coords, :, r))

include("rotations.jl")
include("bonds.jl")