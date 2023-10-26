include("utils.jl")
include("residue_array.jl")

const Backbone = ResidueArray{3}
const BackboneAndOxygen = ResidueArray{4}

nitrogen_coord_matrix(ra::ResidueArray)    = _ith_column_coords(ra, 1)
alphacarbon_coord_matrix(ra::ResidueArray) = _ith_column_coords(ra, 2)
carbon_coord_matrix(ra::ResidueArray)      = _ith_column_coords(ra, 3)
oxygen_coord_matrix(ra::ResidueArray)      = _ith_column_coords(ra, 4)

nitrogens(ra::ResidueArray)    = eachcol(_ith_column_coords(ra, 1))
alphacarbons(ra::ResidueArray) = eachcol(_ith_column_coords(ra, 2))
carbons(ra::ResidueArray)      = eachcol(_ith_column_coords(ra, 3))
oxygens(ra::ResidueArray)      = eachcol(_ith_column_coords(ra, 4))

unroll_coords(ra::ResidueArray{A}) where A = reshape(permutedims(ra.coords, (2, 1, 3)), :, A)

include("oxygen.jl")
include("residue.jl")