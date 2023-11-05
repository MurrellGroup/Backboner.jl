function get_rotation_matrix(
    point1::AbstractVector{T}, center::AbstractVector{T}, point3::AbstractVector{T}
) where T <: Real
    direction1 = point3 - center
    direction2 = point1 - center
    basis1 = normalize(direction1)
    orthogonal_component = direction2 - basis1 * (basis1' * direction2)
    basis2 = normalize(orthogonal_component)
    basis3 = cross(basis1, basis2)
    rotation_matrix = [basis1 basis2 basis3]
    return rotation_matrix
end

# the average rotation matrix calculated from a set of PDB files
const magic_vector = [-0.672, -1.034, 0.003]

function estimate_oxygen_position(
    CA::AbstractVector{T}, C::AbstractVector{T}, next_N::AbstractVector{T},
) where T <: Real
    R = get_rotation_matrix(CA, C, next_N)
    return R' \ magic_vector + C
end

"""
    chain_with_oxygen(chain::Chain{3})

Returns a chain with oxygen atoms added to the residues.
"""
function chain_with_oxygen(
    chain::Chain{3, T},
) where T <: Real
    L = length(chain)

    CAs = eachcol(alphacarbon_coord_matrix(chain))
    Cs = eachcol(carbon_coord_matrix(chain))
    next_Ns = view(eachcol(nitrogen_coord_matrix(chain)), 2:L)

    oxygen_coords = zeros(T, 3, L-1)
    for (i, (CA, C, next_N)) in enumerate(zip(CAs, Cs, next_Ns))
        oxygen_coords[:, i] = estimate_oxygen_position(CA, C, next_N)
    end

    coords = zeros(T, 3, 4, L-1)
    coords[:, 1:3, :] = chain.coords[:, :, 1:L-1]
    coords[:, 4, :] = oxygen_coords

    return Chain(chain.id, coords)
end

"""
    backbone_with_oxygen(backbone::Backbone{3})

Returns a backbone with oxygen atoms added to the residues.
"""
function backbone_with_oxygen(backbone::Backbone{3, T}) where T
    new_chains = Vector{Chain{4,T}}()
    for chain in backbone
        push!(new_chains, chain_with_oxygen(chain))
    end
    return Backbone(new_chains)
end