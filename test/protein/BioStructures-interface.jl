@testset "BioStructures-interface.jl" begin

    @testset "PDB" begin

        @testset "read" begin
            chains_pdb = readpdb("data/1ASS.pdb")
            @test length.(chains_pdb) == [152]
            chains_cif = readchains("data/1ASS.cif", MMCIFFormat)
            @test chains_pdb[1].backbone.coords == chains_cif[1].backbone.coords
        end

        @testset "write" begin
            chains = readpdb("data/1ASS.pdb")
            new_chains = mktempdir() do temp_dir
                temp_path = joinpath(temp_dir, "temp.pdb")
                writepdb(chains, temp_path)
                readpdb(temp_path)
            end
            @test chains[1].backbone.coords == new_chains[1].backbone.coords
        end

    end

end