export Backbone, atom_coord_matrix

"""
    Backbone{A, T}

A wrapper for a 3xAxN matrix of coordinates of atoms in a backbone chain.
"""
struct Backbone{A, T <: Real} <: AbstractArray{T, 3}
    coords::AbstractArray{T, 3}

    function Backbone(coords::AbstractArray{T, 3}) where T <: Real
        A = size(coords, 2)
        @assert size(coords, 1) == 3 "coords must have 3 coordinates per atom"
        return new{A, T}(coords)
    end
end

@inline Base.size(bb::Backbone) = size(bb.coords)
@inline Base.length(bb::Backbone) = size(bb, 3)
@inline Base.getindex(bb::Backbone, i, j, k) = bb.coords[i,j,k]
@inline Base.getindex(bb::Backbone, r::UnitRange{<:Integer}) = Backbone(view(bb.coords, :, :, r))
@inline Base.getindex(bb::Backbone, i::Integer) = view(bb.coords, :, :, i)

"""
    atom_coord_matrix(bb, i)

Returns the coordinates of specific columns of atoms in a backbone.
"""
function atom_coord_matrix(bb::Backbone{A}, i) where A
    return view(bb.coords, :, i, :)
end

"""
    unroll_atoms(bb, i)

Returns the coordinates of specific columns of atoms in a backbone,
but unrolled into a 3xN matrix where N is the number of residues times
the number of columns selected (atoms selected per residue).
"""
function unroll_atoms(bb::Backbone{A}, i=Colon()) where A
    return reshape(atom_coord_matrix(bb, i), 3, :)
end