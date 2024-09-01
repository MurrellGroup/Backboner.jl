module Backboner

export Backbone, coords

export ChainedBonds
export append_bonds!, append_bonds
export prepend_bonds!, prepend_bonds
export get_bond_vectors
export get_bond_lengths
export get_bond_angles
export get_torsion_angles

export Frames

export is_knotted

include("backbone.jl")
include("bonds.jl")
include("frames.jl")
include("knots.jl")

using PrecompileTools

@compile_workload begin
    backbone = Backbone(randn(3, 6))

    bonds = ChainedBonds(backbone)
    Backbone(bonds)

    standard = randn(3,3)
    frames = Frames(backbone, standard)
    Backbone(frames, standard)

    is_knotted(backbone)
end

end