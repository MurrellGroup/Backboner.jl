export oxygen_coords

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
    return R' \ magic_vector + C # inverse of R' * (oxygen_pos - C)
end

# last oxygen is put on the same plane as the last residue triangle
# 1.2 angstroms away from the C atom
# with a 120 degree angle between the C-Ca and C-O bonds
function get_last_oxygen(
    N_pos::V, Ca_pos::V, C_pos::V,
) where V <: AbstractVector{T} where T <: Real
    angle = 2π/3
    bond_length = 1.2
    v = Ca_pos - C_pos
    w = cross(v, Ca_pos - N_pos)
    u = cross(w, v)
    O_pos = C_pos + normalize(u)*bond_length*cos(angle - 0.5π) - normalize(v)*bond_length*sin(angle - 0.5π)

    return O_pos
end

"""
    oxygen_coords(chain::Chain)

Add oxygen atoms to the backbone of a protein, turning the coordinate array from size 3x3xL to 3x4xL-1,
where L is the length of the backbone.
"""
function oxygen_coords(
    chain::Chain,
)
    backbone = chain.backbone
    T = eltype(eltype(backbone))
    L = length(chain)

    CAs = eachcol(alphacarbon_coords(chain))
    Cs = eachcol(carbonyl_coords(chain))
    next_Ns = eachcol(@view(nitrogen_coords(chain)[:, 2:end]))

    oxygen_coords = zeros(T, 3, L)
    for (i, (CA, C, next_N)) in enumerate(zip(CAs, Cs, next_Ns))
        oxygen_coords[:, i] = estimate_oxygen_position(CA, C, next_N)
    end
    oxygen_coords[:, end] = get_last_oxygen(collect(backbone[end-2:end])...)

    return oxygen_coords
end

NCaCO_coords(chain::Chain) = cat(reshape(chain.backbone.coords, 3, 3, :), reshape(oxygen_coords(chain), 3, 1, :), dims=2)