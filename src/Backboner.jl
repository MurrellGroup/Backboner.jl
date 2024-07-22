module Backboner

export Backbone, coords

export ChainedBonds
export append_bonds!, append_bonds
export prepend_bonds!, prepend_bonds
export get_bond_vectors
export get_bond_lengths
export get_bond_angles
export get_dihedrals

export Frames

export is_knotted

export idealize

export Protein

include("backbone.jl")
include("bonds.jl")
include("frames.jl")
include("knots.jl")
include("idealization.jl")
include("protein/protein.jl")

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