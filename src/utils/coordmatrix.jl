export
    atom_coord_matrix,
    unroll_atoms

"""
    atom_coord_matrix(chain, i)

Returns the coordinates of specific columns of atoms in a chain.
"""
function atom_coord_matrix(chain::AbstractChain, i)
    return view(chain.coords, :, i, :)
end

"""
    unroll_atoms(chain, i)

Returns the coordinates of specific columns of atoms in a chain,
but unrolled into a 3xN matrix where N is the number of residues times
the number of columns selected (atoms selected per residue).
"""
function unroll_atoms(chain::AbstractChain, i=:)
    return reshape(atom_coord_matrix(chain, i), 3, :)
end

export
    nitrogen_coord_matrix,
    alphacarbon_coord_matrix,
    carbon_coord_matrix,
    oxygen_coord_matrix

nitrogen_coord_matrix(chain::AbstractChain) = atom_coord_matrix(chain, 1)
alphacarbon_coord_matrix(chain::AbstractChain) = atom_coord_matrix(chain, 2)
carbon_coord_matrix(chain::AbstractChain) = atom_coord_matrix(chain, 3)
oxygen_coord_matrix(chain::AbstractChain) = atom_coord_matrix(chain, 4)