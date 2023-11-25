@testset "residue.jl" begin
    backbone = Backbone(randn(3, 4, 1))
    residue = Residue(1, backbone, 'A', Helix)

    @test summary(residue) == "Residue 1 ALA  Helix"

    io = IOBuffer()
    show(io, residue)
    @test String(take!(io)) == summary(residue)
end