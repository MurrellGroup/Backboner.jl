export
    atom_coord_matrix,
    nitrogen_coord_matrix,
    alphacarbon_coord_matrix,
    carbon_coord_matrix,
    oxygen_coord_matrix

"""
    atom_coord_matrix(backbone, i)

Returns the coordinates of specific columns of atoms in a backbone.
"""
function atom_coord_matrix(backbone::Backbone, i)
    return view(backbone.coords, :, i, :)
end

nitrogen_coord_matrix(backbone::Backbone) = atom_coord_matrix(backbone, 1)
alphacarbon_coord_matrix(backbone::Backbone) = atom_coord_matrix(backbone, 2)
carbon_coord_matrix(backbone::Backbone) = atom_coord_matrix(backbone, 3)
oxygen_coord_matrix(backbone::Backbone) = atom_coord_matrix(backbone, 4)