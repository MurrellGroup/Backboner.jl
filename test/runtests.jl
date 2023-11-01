using BackBoner
using Test

@testset "BackBoner.jl" begin

    @testset "secondary_structure.jl" begin
        @test Loop == SecondaryStructure(1)
        @test Helix == SecondaryStructure(2)
        @test Strand == SecondaryStructure(3)
    end

    include("chains.jl")
    include("backbone.jl")
    include("io.jl")
    include("assign.jl")

    @testset "ncaco.jl" begin
        chain = Chain("A", randn(3, 4, 3))
        @test nitrogen_coord_matrix(chain) == chain.coords[:, 1, :]
        @test alphacarbon_coord_matrix(chain) == chain.coords[:, 2, :]
        @test carbon_coord_matrix(chain) == chain.coords[:, 3, :]
        @test oxygen_coord_matrix(chain) == chain.coords[:, 4, :]
    end

end
