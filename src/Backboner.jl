module Backboner

using LinearAlgebra

include("backbone.jl")
include("bonds.jl")
include("protein/protein.jl")

using .Protein

end