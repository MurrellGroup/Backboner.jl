module BackBoner

export ResidueArray, Backbone, BackboneAndOxygen

using LinearAlgebra

include("AssigningSecondaryStructure/AssigningSecondaryStructure.jl")
import .AssigningSecondaryStructure as ASS
export ASS

include("backbone/backbone.jl")
include("visuals/visuals.jl")

end