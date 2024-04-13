@testset "residue.jl" begin
    residue = Residue(1, Atom[], 'A', 'H')

    io = IOBuffer()
    show(io, residue)
    @test String(take!(io)) == "Residue 1 ALA with 0 atoms (H)"
end