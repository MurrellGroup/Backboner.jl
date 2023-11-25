using Backboner
using Test

@testset "Backboner.jl" begin

    include("secondarystructure.jl")
    include("backbone/backbone.jl")
    include("residue.jl")
    include("chain.jl")
    include("protein.jl")
    include("assign.jl")
    include("io.jl")
    include("coordmatrix.jl")

end
