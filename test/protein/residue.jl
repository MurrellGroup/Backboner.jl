@testset "residue.jl" begin
    residue = Residue(1, 'A', 'H')

    io = IOBuffer()
    show(io, residue)
    @test String(take!(io)) == "Residue 1 ALA H"
end