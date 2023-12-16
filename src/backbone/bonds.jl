export atom_distances, carbonyl_nitrogen_distances

function atom_displacements(backbone::Backbone{A}, atom1::Integer, atom2::Integer, residue_offset::Integer=0) where A
    @assert 1 <= atom1 <= A && 1 <= atom2 <= A "Backbone{$N} does not have atoms $atom1 and $atom2"
    coords = backbone.coords
    displacements = @view(coords[:, atom2, 1+residue_offset:end]) .- @view(coords[:, atom1, 1:end-residue_offset])
    return displacements
end

function backbone_bond_vectors(backbone::Backbone{A}) where A
    backbone1 = Backbone{1}(backbone)
    bond_vectors = atom_displacements(backbone1, 1, 1, 1)
    return bond_vectors
end

"""
    atom_distances(backbone::Backbone, atom1::Integer, atom2::Integer, residue_offset::Integer)

Calculate the distances between all pairs of two types atoms in a backbone, e.g.
the distances between all pairs of contiguous carbonyl and nitrogen atoms.
atom1 and atom2 are the indices of the atoms in the backbone, and residue_offset
is the number of residues between the atoms (0 by default).

Returns a vector of distances of length (length(backbone) - residue_offset).
"""
function atom_distances(backbone::Backbone{A}, atom1::Integer, atom2::Integer, residue_offset::Integer=0) where A
    @assert 1 <= atom1 <= A && 1 <= atom2 <= A "Backbone{$N} does not have atoms $atom1 and $atom2"
    displacements = atom_displacements(backbone, atom1, atom2, residue_offset)
    distances = reshape(mapslices(norm, displacements, dims=1), :)
    return distances
end