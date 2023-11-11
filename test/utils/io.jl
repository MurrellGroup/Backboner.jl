@testset "io.jl" begin

    @testset "PDB" begin

        @testset "read" begin
            protein = pdb_to_protein("data/1ASS.pdb")
            @test length.(protein) == [152]
        end

        @testset "write" begin
            protein = pdb_to_protein("data/1ASS.pdb")
            protein_to_pdb(protein, "temp.pdb")
            protein2 = pdb_to_protein("temp.pdb")
            @test protein == protein2
            rm("temp.pdb")
        end

    end

end