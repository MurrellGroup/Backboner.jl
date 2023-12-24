using Backboner
using Test

using LinearAlgebra

@testset "Backboner.jl" begin

    include("backbone.jl")
    include("bonds.jl")
    include("protein/protein.jl")

end
