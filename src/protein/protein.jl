module Protein

# TODO: un-export Atom and Residue as they shouldn't really be public (breaking)

using ..Backboner

import BioStructures

include("atom.jl")
export Atom

include("residue.jl")
export Residue

include("chain.jl")
export Chain
export has_assigned_ss

include("oxygen.jl")
export oxygen_coords

include("BioStructures-interface.jl")
export nitrogen_coords, alphacarbon_coords, carbonyl_coords
export nitrogen_alphacarbon_distances, alphacarbon_carbonyl_distances, carbonyl_nitrogen_distances
export phi_angles, psi_angles, omega_angles

include("idealization.jl")
export readchains, writechains
export readpdb, writepdb
export readcif, writecif
export PDBFormat, MMCIFFormat

PDBEntry(pdbid::AbstractString; format=BioStructures.PDBFormat, kwargs...) = mktempdir() do dir
    path = BioStructures.downloadpdb(pdbid; dir, format, kwargs...)
    Backboner.Protein.readchains(path, format)
end

macro pdb_str(pdbid)
    quote
        PDBEntry($(esc(pdbid)))
    end
end

export PDBEntry, @pdb_str

end