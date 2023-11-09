using Backboner
using Test

@testset "Backboner.jl" begin

    @testset "secondary_structure.jl" begin
        @test Loop == SecondaryStructure(1)
        @test Helix == SecondaryStructure(2)
        @test Strand == SecondaryStructure(3)
    end

    include("chains/chains.jl")
    include("utils/utils.jl")
    include("backbone.jl")
    include("assign.jl")
    include("io.jl")

end
