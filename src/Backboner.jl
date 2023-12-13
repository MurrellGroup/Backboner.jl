module Backboner

using LinearAlgebra
using Flux:batched_mul

import Rotations
import PDBTools

export has_assigned_ss

has_assigned_ss(ssvector::Vector{Char}) = all(!=(' '), ssvector)

include("backbone/backbone.jl")
include("residue.jl")
include("chain.jl")
include("protein.jl")
include("assign.jl")
include("io.jl")
include("dihedrals.jl")

end