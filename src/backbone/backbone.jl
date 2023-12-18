export Backbone, atom_coords

"""
    Backbone{A, T <: Real} <: AbstractArray{T, 3}

The `Backbone{A}` type is designed to efficiently store and manipulate the three-dimensional coordinates of backbone atoms.
The type parameter `A` stands for the number of atoms per residue, and is used in the size of the coordinate array, 3xAxL,
where L is the number of residues.

Backbone{3} is used to store 3-dimensional coordinates of the contiguous backbone atoms (N, Cα, C) of a protein chain.
"""
struct Backbone{A, T <: Real} <: AbstractArray{T, 3}
    coords::AbstractArray{T, 3}

    function Backbone{A}(coords::AbstractArray{T, 3}) where {A, T}
        @assert size(coords, 1) == 3 "coords must have 3 coordinates per atom"
        @assert size(coords, 2) == A "The second dimension of coords must be $A"
        return new{A, T}(coords)
    end

    Backbone(coords::AbstractArray{T, 3}) where T = Backbone{size(coords, 2)}(coords)

    Backbone{A}(backbone::Backbone) where A = Backbone{A}(reshape(backbone.coords, 3, A, :))
end

@inline Base.size(backbone::Backbone) = size(backbone.coords)
@inline Base.length(backbone::Backbone) = size(backbone, 2) * size(backbone, 3)
@inline Base.getindex(backbone::Backbone, i, j, k) = backbone.coords[i, j, k]
@inline Base.getindex(backbone::Backbone, i::Integer) = view(backbone.coords, :, :, i)
@inline Base.getindex(backbone::Backbone, r::UnitRange{Int}) = Backbone(view(backbone.coords, :, :, r))

"""
    atom_coords(backbone, i)

Returns the coordinates of specific columns of atoms in a backbone.
"""
@inline atom_coords(backbone::Backbone, i) = view(backbone.coords, :, i, :)

include("rotations.jl")
include("bonds.jl")