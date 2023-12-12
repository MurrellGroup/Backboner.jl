export atom_distances, carbonyl_nitrogen_distances

function atom_distances(backbone::Backbone{N}, atom1::Integer, atom2::Integer, residue_offset::Integer) where N
    @assert 1 <= atom1 <= N && 1 <= atom2 <= N "Backbone{$N} does not have atoms $atom1 and $atom2"
    coords = backbone.coords
    displacements = coords[:, atom1, 1:end-residue_offset] .- coords[:, atom2, 1+residue_offset:end]
    distances = dropdims(mapslices(norm, displacements, dims=1), dims=1)
    return distances
end

carbonyl_nitrogen_distances(backbone::Backbone) = atom_distances(backbone, 3, 1, 1)