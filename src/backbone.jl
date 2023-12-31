export Backbone

"""
    Backbone{T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{AbstractVector{T}}

The `Backbone` type is designed to efficiently store and manipulate the three-dimensional coordinates of backbone atoms.

# Examples

A `Backbone` can be created from a matrix of coordinates:

```jldoctest
julia> backbone = Backbone(zeros(3, 5)) # 5 atoms with 3 coordinates each
5-element Backbone{Float64}:
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]

julia> backbone[1] = [1.0, 2.0, 3.0]; # set the first atom's coordinates

julia> backbone
5-element Backbone{Float64}:
 [1.0, 2.0, 3.0]
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]

julia> backbone[1:2] # indexing by range returns a `Backbone` wrapped around a view of the original
2-element Backbone{Float64}:
 [1.0, 2.0, 3.0]
 [0.0, 0.0, 0.0]
```

Arrays will always be flattened to be a 3xN matrix:

```jldoctest
julia> backbone = Backbone(zeros(3, 3, 100)) # 3 coordinates per atom, 3 atoms per residue, 100 residues
300-element Backbone{Float64}:
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]
 ⋮
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]
```
"""
struct Backbone{T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{AbstractVector{T}}
    coords::M

    function Backbone{T, M}(coords::M) where {T <: Real, M <: AbstractMatrix{T}}
        @assert size(coords, 1) == 3 "coords must have 3 coordinates per atom"
        return new{T, M}(coords)
    end

    function Backbone{T}(coords::AbstractArray{T}) where T
        new_coords = reshape(coords, size(coords, 1), :)
        M = typeof(new_coords)
        return Backbone{T, M}(new_coords)
    end
end

Backbone(coords::AbstractArray{T}) where T = Backbone{T}(coords)

Backbone{T}(::UndefInitializer, n_atoms::Integer) where T = Backbone{T}(Matrix{T}(undef, 3, n_atoms))

@inline Base.:(==)(backbone1::Backbone, backbone2::Backbone) = backbone1.coords == backbone2.coords
@inline Base.:(≈)(backbone1::Backbone, backbone2::Backbone) = backbone1.coords ≈ backbone2.coords

@inline Base.length(backbone::Backbone) = size(backbone.coords, 2)
@inline Base.size(backbone::Backbone) = Tuple(length(backbone))

@inline Base.getindex(backbone::Backbone, i::Integer) = backbone.coords[:, i]
@inline Base.getindex(backbone::Backbone, r::AbstractVector{<:Integer}) = Backbone(backbone.coords[:, r])
@inline Base.view(backbone::Backbone, I...) = Backbone(view(backbone.coords, :, I...))

@inline Base.setindex!(backbone::Backbone, coords::AbstractVector, i::Integer) = (backbone.coords[:, i] .= coords)
@inline Base.setindex!(backbone::Backbone, coords::AbstractMatrix, r::AbstractVector{<:Integer}) = (backbone.coords[:, r] .= coords)
