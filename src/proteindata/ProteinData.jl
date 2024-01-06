module ProteinData

using ..Backboner

nitrogen_alphacarbon_distances(backbone::Backbone) = get_atom_distances(backbone, 1, 1, 3)
alphacarbon_carbonyl_distances(backbone::Backbone) = get_atom_distances(backbone, 2, 1, 3)
carbonyl_nitrogen_distances(backbone::Backbone) = get_atom_distances(backbone, 3, 1, 3)

include("numbering.jl")
include("pdb.jl")
include("records.jl")

export GeneralRecord

end