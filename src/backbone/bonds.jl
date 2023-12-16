export get_atom_displacements, get_atom_distances, get_bond_vectors, get_bond_angles

function get_atom_displacements(backbone::Backbone{A}, atom1::Integer, atom2::Integer, residue_offset::Integer=0) where A
    @assert 1 <= atom1 <= A && 1 <= atom2 <= A "Backbone{$N} does not have atoms $atom1 and $atom2"
    coords = backbone.coords
    displacements = @view(coords[:, atom2, 1+residue_offset:end]) .- @view(coords[:, atom1, 1:end-residue_offset])
    return displacements
end

norms(vectors::AbstractMatrix) = reshape(mapslices(norm, vectors, dims=1), :)

"""
    get_atom_distances(backbone::Backbone, atom1::Integer, atom2::Integer, residue_offset::Integer)

Calculate the distances between all pairs of two types atoms in a backbone, e.g.
the distances between all pairs of contiguous carbonyl and nitrogen atoms.
atom1 and atom2 are the indices of the atoms in the backbone, and residue_offset
is the number of residues between the atoms (0 by default).

Returns a vector of distances of length (length(backbone) - residue_offset).
"""
function get_atom_distances(backbone::Backbone{A}, atom1::Integer, atom2::Integer, residue_offset::Integer=0) where A
    @assert 1 <= atom1 <= A && 1 <= atom2 <= A "Backbone{$N} does not have atoms $atom1 and $atom2"
    displacements = get_atom_displacements(backbone, atom1, atom2, residue_offset)
    distances = norms(displacements)
    return distances
end


function get_bond_vectors(backbone::Backbone{A}) where A
    backbone1 = Backbone{1}(backbone)
    bond_vectors = get_atom_displacements(backbone1, 1, 1, 1)
    return bond_vectors
end


function get_bond_lengths(bond_vectors::AbstractVector{<:AbstractVector{T}}) where T
    bond_lengths = norms(bond_vectors)
    return bond_lengths
end

get_bond_lengths(backbone::Backbone) = get_bond_lengths(get_bond_vectors(backbone))


get_bond_angle(v1::T, v2::T) where T <: AbstractVector = acos(dot(v1, v2) / (norm(v1) * norm(v2)))

function get_bond_angles(bond_vectors::AbstractMatrix{T}) where T
    bond_vector_pairs = zip(eachcol(bond_vectors), Iterators.drop(eachcol(bond_vectors), 1))
    bond_angles = [get_bond_angle(v1, v2) for (v1, v2) in bond_vector_pairs]
    return bond_angles
end

get_bond_angles(backbone::Backbone) = get_bond_angles(get_bond_vectors(backbone))


struct ContinuousBonds{T} <: AbstractVector{AbstractVector{T}}
    vectors::AbstractMatrix{T}
    lengths::AbstractVector{T}
    angles::AbstractVector{T}

    function ContinuousBonds(vectors::AbstractMatrix{T}) where T
        lengths = get_bond_lengths(vectors)
        angles = get_bond_angles(vectors)

        L = size(vectors, 2) + 1
        @assert size(lengths) == L-1
        @assert size(angles) == L-2

        return new{T}(vectors, lengths, angles)
    end

    function ContinuousBonds(backbone::Backbone{A, T}) where {A, T}
        return ContinuousBonds(get_bond_vectors(backbone))
    end
end

Base.length(bonds::ContinuousBonds) = size(bonds.vectors, 2)
Base.size(bonds::ContinuousBonds) = Tuple(length(bonds))
Base.getindex(bonds::ContinuousBonds, i) = @view(bonds.vectors[:, i])