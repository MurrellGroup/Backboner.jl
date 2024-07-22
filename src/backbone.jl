"""
    Backbone{T<:Real,M<:AbstractMatrix{T}} <: AbstractVector{AbstractVector{T}}

The `Backbone` type is designed to efficiently store and manipulate the three-dimensional coordinates of backbone atoms.
"""
struct Backbone{T<:Real,M<:AbstractMatrix{T}} <: AbstractVector{AbstractVector{T}}
    coords::M

    function Backbone{T,M}(coords::M) where {T,M}
        size(coords, 1) == 3 || throw(ArgumentError("Expected the first dimension of coords to have a size of 3"))
        return new(coords)
    end
end

Backbone(coords::M) where {T<:Real,M<:AbstractMatrix{T}} = Backbone{T,M}(coords)
Backbone(coords::A) where {T<:Real,A<:AbstractArray{T}} = Backbone(reshape(coords, size(coords, 1), :))

coords(backbone::Backbone) = backbone.coords

@inline Base.size(backbone::Backbone) = Tuple(size(backbone.coords, 2))
@inline Base.getindex(backbone::Backbone, i) = Backbone(backbone.coords[:, i])
@inline Base.view(backbone::Backbone, i) = Backbone(view(backbone.coords, :, i))
@inline Base.setindex!(backbone::Backbone, v, i) = (backbone.coords[:, i] .= v)

@inline Base.getindex(backbone::Backbone, i::Integer) = backbone.coords[:, i]
@inline Base.view(backbone::Backbone, i::Integer) = view(backbone.coords, :, i)

Base.:(==)(backbone1::Backbone, backbone2::Backbone) = backbone1.coords == backbone2.coords

Base.hash(backbone::Backbone, h::UInt) = hash(backbone.coords, h)