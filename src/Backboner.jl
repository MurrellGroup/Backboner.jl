module Backboner

using LinearAlgebra
using Rotations: AngleAxis
using NNlib: batched_mul

norms(A::AbstractArray{<:Number}; dims=1) = sqrt.(sum(abs2, A; dims))
dots(A1::AbstractArray{<:Number}, A2::AbstractMatrix{<:Number}; dims=1) = sum(A1 .* A2; dims)
normalize_slices(A::AbstractArray{<:Number}; dims=1) = A ./ norms(A; dims)

include("backbone.jl")
export Backbone, coords

include("bonds.jl")
export ChainedBonds
export append_bonds!, append_bonds
export prepend_bonds!, prepend_bonds
export get_bond_vectors
export get_bond_lengths
export get_bond_angles
export get_torsion_angles

include("frames.jl")
export Frames

include("knots.jl")
export is_knotted

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