module Protein

# TODO: un-export Atom and Residue as they shouldn't really be public (breaking)

export Atom

export Residue

export Chain
export has_assigned_ss

export nitrogen_coords, alphacarbon_coords, carbonyl_coords
export nitrogen_alphacarbon_distances, alphacarbon_carbonyl_distances, carbonyl_nitrogen_distances
export phi_angles, psi_angles, omega_angles

export readchains, writechains
export readpdb, writepdb
export readcif, writecif
export PDBFormat, MMCIFFormat

using ..Backboner

import BioStructures

include("atom.jl")
include("residue.jl")
include("chain.jl")
include("oxygen.jl")
include("BioStructures-interface.jl")
include("idealization.jl")

end