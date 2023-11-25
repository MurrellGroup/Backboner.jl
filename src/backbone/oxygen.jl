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

function add_last_oxygen!(
    coords::AbstractArray{T, 3},
) where T
    @assert size(coords, 1) == 3
    @assert size(coords, 2) == 4
    @assert size(coords, 3) > 1

    L = size(coords, 3)

    prev_N, prev_CA, prev_C, prev_O = eachcol(coords[:, :, L-1])
    last_N, last_CA, last_C, last_O = eachcol(coords[:, :, L])

    prev_rot_matrix = get_rotation_matrix(prev_N, prev_CA, prev_C)
    last_rot_matrix = get_rotation_matrix(last_N, last_CA, last_C)

    next_N_relative_to_triangle = prev_rot_matrix' * (last_N - prev_CA)

    OXT = last_rot_matrix' \ next_N_relative_to_triangle + last_CA
    last_O = estimate_oxygen_position(last_CA, last_C, OXT)
    coords[:, 4, L] = last_O

    return coords
end

"""
    add_oxygens(backbone::Backbone{3})

Add oxygen atoms to the backbone of a protein, turning the coordinate array from size 3x3xL to 3x4xL-1,
where L is the length of the backbone.

!!! note 
    One residue is lost in the process, since the orientation of the last oxygen atom cannot be determined.
    We may consider adding a feature for creating a dummy oxygen atom at the end of the backbone with
    randomized orientation to preserve the length of the backbone.
"""
function add_oxygens(
    backbone::Backbone{3,T},
) where T <: Real
    L = length(backbone)

    CAs = eachcol(alphacarbon_coord_matrix(backbone))
    Cs = eachcol(carbon_coord_matrix(backbone))
    next_Ns = eachcol(nitrogen_coord_matrix(backbone[2:end]))

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