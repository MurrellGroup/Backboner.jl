export Backbone, remove_column

"""
    Backbone{N, T <: Real} <: AbstractArray{T, 3}

A wrapper for a 3xNxL array of coordinates of atoms.
Backbone{4} is used to store 3-dimensional coordinates of the backbone atoms (N, CA, C, O) of a protein chain.
"""
struct Backbone{N, T <: Real} <: AbstractArray{T, 3}
    coords::AbstractArray{T, 3}

    function Backbone{N}(coords::AbstractArray{T, 3}) where {N, T}
        @assert size(coords, 1) == 3 "coords must have 3 coordinates per atom"
        @assert size(coords, 2) == N "The second dimension of coords must be $N"
        return new{N, T}(coords)
    end

    Backbone(coords::AbstractArray{T, 3}) where T = Backbone{size(coords, 2)}(coords)
end

@inline Base.:(==)(backbone1::Backbone, backbone2::Backbone) = backbone1.coords == backbone2.coords
@inline Base.size(backbone::Backbone) = size(backbone.coords)
@inline Base.length(backbone::Backbone) = size(backbone, 3)
@inline Base.getindex(backbone::Backbone, i, j, k) = backbone.coords[i,j,k]
@inline Base.getindex(backbone::Backbone, i::Integer) = view(backbone.coords, :, :, i)
@inline Base.getindex(backbone::Backbone, r::UnitRange{Int}) = Backbone(view(backbone.coords, :, :, r))

function remove_column(backbone::Backbone{N, T}, i::Integer) where {N, T}
    @assert 1 <= i <= N "i must be between 1 and $N"
    Backbone(view(backbone.coords, :, [1:i-1;i+1:N], :))
end

"""
    atom_coord_matrix(backbone, i)

Returns the coordinates of specific columns of atoms in a backbone.
"""
function atom_coord_matrix(backbone::Backbone, i)
    return view(backbone.coords, :, i, :)
end

include("rotations.jl")
include("oxygen.jl")