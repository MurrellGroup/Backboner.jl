using Backboner
using Test

using LinearAlgebra

@testset "Backboner.jl" begin

    include("backbone.jl")
    include("bonds.jl")
    include("frames.jl")
    include("knots.jl")
    include("protein/protein.jl")

end
