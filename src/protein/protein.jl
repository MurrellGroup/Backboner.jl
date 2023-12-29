module Protein

using ..Backboner
using LinearAlgebra

import PDBTools

include("residue.jl")
include("chain.jl")
include("oxygen.jl")
include("assign.jl")
include("io.jl")
include("rotations.jl")

end