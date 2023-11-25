@testset "residue.jl" begin
    backbone = Backbone(randn(3, 4, 1))
    residue = Residue(1, backbone, 'A', Helix)

    # index field get padded to digits in length of backbone, ss gets padded to 6 cause Strand is 6 chars long.
    # Unassigned is longer but that's fine.
    @test summary(residue) == "Residue 1 ALA  Helix"

    io = IOBuffer()
    show(io, residue)
    @test String(take!(io)) == summary(residue)
end