@testset "residue.jl" begin
    backbone = Backbone(randn(3, 4, 100))
    residue = Residue(1, backbone, 'A', 'H')
    @test backbone_atom_coords(residue) == backbone[1]

    # index field get padded to digits in length of backbone
    @test summary(residue) == "Residue   1 ALA H"

    io = IOBuffer()
    show(io, residue)
    @test String(take!(io)) == summary(residue)
end