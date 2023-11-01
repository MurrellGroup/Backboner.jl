export nitrogen_coord_matrix, alphacarbon_coord_matrix, carbon_coord_matrix, oxygen_coord_matrix

nitrogen_coord_matrix(chain::Chain) = atom_coord_matrix(chain, 1)
alphacarbon_coord_matrix(chain::Chain) = atom_coord_matrix(chain, 2)
carbon_coord_matrix(chain::Chain) = atom_coord_matrix(chain, 3)
oxygen_coord_matrix(chain::Chain) = atom_coord_matrix(chain, 4)