module Protein

using ..Backboner
import BioStructures

include("residue.jl")
include("chain.jl")
include("oxygen.jl")
include("pdb.jl")

end