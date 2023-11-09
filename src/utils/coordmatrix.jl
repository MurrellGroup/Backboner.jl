export
    atom_coord_matrix,
    nitrogen_coord_matrix,
    alphacarbon_coord_matrix,
    carbon_coord_matrix,
    oxygen_coord_matrix

"""
    atom_coord_matrix(chain, i)

Returns the coordinates of specific columns of atoms in a chain.
"""
function atom_coord_matrix(chain::AbstractChain, i)
    return view(chain.coords, :, i, :)
end

nitrogen_coord_matrix(chain::AbstractChain) = atom_coord_matrix(chain, 1)
alphacarbon_coord_matrix(chain::AbstractChain) = atom_coord_matrix(chain, 2)
carbon_coord_matrix(chain::AbstractChain) = atom_coord_matrix(chain, 3)
oxygen_coord_matrix(chain::AbstractChain) = atom_coord_matrix(chain, 4)