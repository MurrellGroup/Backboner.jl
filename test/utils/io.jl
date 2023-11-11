@testset "io.jl" begin

    @testset "PDB" begin

        @testset "read" begin
            protein = pdb_to_protein("data/1ASS.pdb")
            @test length.(protein) == [152]
        end

        @testset "write" begin
            out = "temp.pdb"
            try
                protein = pdb_to_protein("data/1ASS.pdb")
                protein_to_pdb(protein, out)
                protein2 = pdb_to_protein(out)
                @test protein == protein2
            finally
                rm(out, force=true)
            end
        end

    end

end