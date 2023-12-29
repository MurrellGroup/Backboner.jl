# This algorithm is presented in a paper by Cui et al.: https://almob.biomedcentral.com/articles/10.1186/1748-7188-8-5

import Distances: euclidean, sqeuclidean
import Statistics: mean

# written as D(P_i, P_j) in the paper
bottleneck_distance(P_i::Backbone, P_j::Backbone) = sqrt(maximum(sqeuclidean, zip(P_i, P_j)))
 
# the root-mean-square deviation of atomic positions: https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
RMSD(m1::M, m2::M) where M <: AbstractMatrix{<:Real} = sqrt(mean(_pairwise_column_distances(m1, m2, sqeuclidean)))
RMSD(v1::V, v2::V) where V <: AbstractVector{<:Real} = sqrt(mean(x -> x^2, v1 .- v2)) # should be same as reshaping into matrix

# written as S_f(P_i) in the paper
function free_energy_score(P_i::Backbone) # pass Chain instead (or amino acid vector)
    return 0 # TODO
end

function D_α(P_i::Backbone, P_j::Backbone)
    return RMSD(alphacarbon_coords(P_i), alphacarbon_coords(P_j))
end

angle_D(a, b) = min(abs(a - b), 2π - abs(a - b))

function D_φψ(P_i::Backbone, P_j::Backbone)
    bonds_i = ChainedBonds(P_i)
    bonds_j = ChainedBonds(P_j)
    φ_i, ψ_i = phi_angles(bonds_i), psi_angles(bonds_i)
    φ_j, ψ_j = phi_angles(bonds_j), psi_angles(bonds_j)
    return sqrt(mean(angle_D.(φ_i, φ_j).^2 .+ angle_D.(ψ_i, ψ_j).^2))
end

function D_O(P_i::Backbone, P_j::Backbone)
    return sqrt(RMSD(oxygen_coords(P_i)^2, oxygen_coords(P_j))^2)
end

# TODO: add beta and hydrogen terms
# TODO: precomputed bonds_0
# written as S_{BB}(P_i) in the paper
# evaluates not only the similarity between idealized and target structure,
# but also evaluates the free energy to ensure it's protein-like
function score_function(P_0::Backbone, P_i::Backbone, w1=1, w2=1, w3=1)
    - w1*D_α(P_0, P_i) +
    - w2*D_O(P_0, P_i) +
    - w3*D_φψ(P_0, P_i)
end

# Main Dynamic Programming Algorithm
function idealize_protein_structure(P_0::Backbone, r=1, ϵ=r/5)
    
end




"""
    idealize(backbone, idealized_lengths, idealized_angles)

`idealized_lengths` and `idealized_angles` will be *repeated* indefinitely to fill the backbone.
Thus, for protein chains, 3 elements will suffice for each vector.
"""
function idealize(
    backbone::Backbone{T},
    idealized_lengths::AbstractVector{T},
    idealized_angles::AbstractVector{T},
)

end