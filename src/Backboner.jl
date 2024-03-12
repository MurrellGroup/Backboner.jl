module Backboner

include("backbone.jl")
include("bonds.jl")
include("frames.jl")
include("knots.jl")
include("idealization.jl")
include("protein/protein.jl")

export Protein

using PrecompileTools

@time @compile_workload begin
    backbone = Backbone{Float64}(reshape(1:18, 3, :))

    bonds = ChainedBonds(backbone)
    Backbone(bonds)

    frames = Frames(backbone, Protein.STANDARD_RESIDUE_ANGSTROM)
    Backbone(frames, Protein.STANDARD_RESIDUE_ANGSTROM)

    is_knotted(backbone)
end

end