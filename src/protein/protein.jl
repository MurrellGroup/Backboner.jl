module Protein

using ..Backboner
import BioStructures

include("atom.jl")
include("residue.jl")
include("chain.jl")
include("oxygen.jl")
include("pdb.jl")
include("idealization.jl")

end