module Backboner

include("backbone.jl")
include("bonds.jl")
include("frames.jl")
include("knots.jl")
include("idealization.jl")
include("protein/protein.jl")

export Protein

using PrecompileTools

@compile_workload begin
    backbone = Backbone(randn(3, 6))

    bonds = ChainedBonds(backbone)
    Backbone(bonds)

    frames = Frames(backbone, Protein.STANDARD_TRIANGLE_ANGSTROM)
    Backbone(frames, Protein.STANDARD_TRIANGLE_ANGSTROM)

    is_knotted(backbone)
end

end