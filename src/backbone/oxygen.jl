export add_oxygens

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
function add_last_oxygen!(
    coords::AbstractArray{T, 3},
) where T
    @assert size(coords, 1) == 3
    @assert size(coords, 2) == 4
    @assert size(coords, 3) > 0

    L = size(coords, 3)

    N_pos, Ca_pos, C_pos, _ = eachcol(coords[:, :, L])

    angle = 2π/3
    bond_length = 1.2
    v = Ca_pos - C_pos
    w = cross(v, Ca_pos - N_pos)
    u = cross(w, v)
    O_pos = C_pos + normalize(u)*bond_length*cos(angle - 0.5π) - normalize(v)*bond_length*sin(angle - 0.5π)
    coords[:, 4, L] = O_pos

    return coords
end

"""
    add_oxygens(backbone::Backbone{3})

Add oxygen atoms to the backbone of a protein, turning the coordinate array from size 3x3xL to 3x4xL-1,
where L is the length of the backbone.
"""
function add_oxygens(
    backbone::Backbone{3, T},
) where T <: Real
    L = length(backbone)

    CAs = eachcol(atom_coord_matrix(backbone, 2))
    Cs = eachcol(atom_coord_matrix(backbone, 3))
    next_Ns = eachcol(atom_coord_matrix(backbone[2:end], 1))

    oxygen_coords = zeros(T, 3, L-1)
    for (i, (CA, C, next_N)) in enumerate(zip(CAs, Cs, next_Ns))
        oxygen_coords[:, i] = estimate_oxygen_position(CA, C, next_N)
    end

    coords = zeros(T, 3, 4, L)
    coords[:, 1:3, 1:L] = backbone.coords[:, :, :]
    coords[:, 4, 1:L-1] = oxygen_coords

    add_last_oxygen!(coords)

    return Backbone(coords)
end