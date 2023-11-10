module Backboner

using LinearAlgebra
using Rotations

include("secondarystructure.jl")
include("backbone/backbone.jl")
include("chain/chain.jl")
include("protein.jl")
include("assign.jl")
include("utils/utils.jl")

end