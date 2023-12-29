@testset "residue.jl" begin
    residue = Residue(1, 'A', 'H')

    @test summary(residue) == "Residue 1 ALA H"

    io = IOBuffer()
    show(io, residue)
    @test String(take!(io)) == summary(residue)
end