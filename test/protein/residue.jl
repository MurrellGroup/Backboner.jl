@testset "residue.jl" begin
    residue = Protein.Residue(1, Protein.Atom[], 'A', 'H', 'A')

    io = IOBuffer()
    show(io, residue)
    @test String(take!(io)) == "Residue 1A ALA with 0 atoms (H)"
end