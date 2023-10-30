module BackBoner

using LinearAlgebra

import AssigningSecondaryStructure as ASS
export ASS

include("backbone.jl")
include("oxygen.jl")
include("ncaco.jl")

end