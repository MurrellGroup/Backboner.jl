using Backboner
using Test

@testset "Backboner.jl" begin

    include("secondarystructure.jl")
    include("backbone/backbone.jl")
    include("chain/chain.jl")
    include("protein.jl")
    include("assign.jl")
    include("utils/utils.jl")

end
