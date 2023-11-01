module BackBoner

using LinearAlgebra

include("secondary_structure.jl")
include("chains/chains.jl")
include("backbone.jl")
include("oxygen.jl")
include("ncaco.jl")
include("io.jl")
include("assign.jl")

end