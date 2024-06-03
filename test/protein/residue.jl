@testset "residue.jl" begin
    residue = Residue(1, Atom[], 'A', 'H', 'A')

    io = IOBuffer()
    show(io, residue)
    @test String(take!(io)) == "Residue 1A ALA with 0 atoms (H)"
end