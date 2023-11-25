@testset "residue.jl" begin
    backbone = Backbone(randn(3, 4, 1))
    residue = Residue(1, backbone, 'A', Helix)

    @test summary(chain) == "Residue 1 ALA  Helix"

    io = IOBuffer()
    show(io, chain)
    @test String(take!(io)) == summary(chain)
end