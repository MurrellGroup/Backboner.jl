module Backboner

using LinearAlgebra
using Rotations

import PDBTools

include("secondarystructure.jl")
include("backbone/backbone.jl")
include("chain.jl")
include("protein.jl")
include("assign.jl")
include("io.jl")
include("coordmatrix.jl")

end