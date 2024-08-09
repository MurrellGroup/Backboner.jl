module ZygoteIdealizationExt

using Backboner

import Zygote: Params, withgradient

col_lengths(m::Matrix{<:Real}) = vec(sqrt.(sum(abs2.(m), dims=1)))

function length_loss(
    m::Matrix{T}, ideal::Vector{T};
    mask::Vector{Bool} = ones(Bool, size(m, 2)-1),
) where T <: Real
    lengths = col_lengths(m[:, 2:end] .- m[:, 1:end-1])
    return sum(((lengths .- ideal).^2) .* mask)
end

function col_angles(
    ab::Matrix{T}, cb::Matrix{T},
) where T <: Real
    dots = sum(ab .* cb, dims=1)
    ab_norms = sqrt.(sum(ab.^2, dims=1))
    cb_norms = sqrt.(sum(cb.^2, dims=1))
    angles = acos.(clamp.(dots ./ (ab_norms .* cb_norms), T(-1.0), T(1.0)))
    return vec(angles)
end

function angle_loss(
    m::Matrix{T}, ideal::Vector{T};
    mask::Vector{Bool} = ones(Bool, size(a, 2)-2),
) where T <: Real
    a = m[:, 1:end-2]
    b = m[:, 2:end-1]
    c = m[:, 3:end]
    angles = col_angles(a .- b, c .- b)
    return sum(((sin.(angles) .- sin.(ideal)).^2 .+ (cos.(angles) .- cos.(ideal)).^2) .* mask)
end

function ideal_loss(
    coords::Matrix{T},
    offsets::Matrix{T},
    lengths::Vector{T},
    angles::Vector{T},
    length_mask::Vector{Bool} = ones(Bool, size(original_xyz, 2)-1),
) where T <: Real
    angle_mask = Vector(length_mask[2:end] .& length_mask[1:end-1])
    offset_coords = coords .+ offsets
    offset_loss = sum(abs2.(offsets))
    l_loss = length_loss(offset_coords, lengths, mask = length_mask)
    a_loss = angle_loss(offset_coords, angles, mask = angle_mask)
    return offset_loss, l_loss, a_loss
end


"""
    idealize(
        backbone::Backbone{T},
        ideal_lengths::Vector{T},
        ideal_angles::Vector{T};
        mask_tolerance = 0.5,
        n_iterations = 300,
    ) where T <: Real

Idealizes a `Backbone` by adjusting the coordinates to match the ideal lengths and angles.
"""
function Backboner.idealize(
    backbone::Backbone{T},
    ideal_lengths::Vector{T},
    ideal_angles::Vector{T};
    mask_tolerance = 0.5,
    n_iterations = 300,
) where T <: Real
    bonds = ChainedBonds(backbone)
    ideal_lengths = [ideal for (_, ideal) in zip(get_bond_lengths(bonds), Iterators.cycle(ideal_lengths))]
    ideal_angles = [ideal for (_, ideal) in zip(get_bond_angles(bonds), Iterators.cycle(ideal_angles))]
    length_mask = Vector(abs.(get_bond_lengths(bonds) .- ideal_lengths) .< mask_tolerance)

    offsets = zeros(T, size(backbone.coords))

    function total_loss(coords, offsets, w1, w2, w3)
        offset_loss, l_loss, a_loss = ideal_loss(coords, offsets, ideal_lengths, ideal_angles, length_mask)
        return w1*offset_loss^3 + w2*l_loss + w3*a_loss
    end

    for i in 1:n_iterations
        l, g = withgradient(Params([offsets])) do
            total_loss(backbone.coords, offsets, 0, 1.4, 1) * 0.13
        end

        offsets .-= g.grads[offsets]
    end

    return Backbone(backbone.coords .+ offsets)
end

end