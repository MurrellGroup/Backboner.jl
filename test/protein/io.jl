@testset "io.jl" begin

    @testset "PDB" begin

        @testset "read" begin
            protein = readpdb("data/1ASS.pdb")
            @test length.(protein) == [152]
        end

        @testset "write" begin
            out = "temp.pdb"
            try
                protein = readpdb("data/1ASS.pdb")
                writepdb(protein, out)
                protein2 = readpdb(out)
                @test protein == protein2
            finally
                rm(out, force=true)
            end
        end

    end

end