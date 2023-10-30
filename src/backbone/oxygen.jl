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

const magic_vector = [-0.672, -1.034, 0.003]

function get_oxygen_coords(
    CA::AbstractVector{T}, C::AbstractVector{T}, next_N::AbstractVector{T},
) where T <: Real
    R = get_rotation_matrix(CA, C, next_N)
    return R' \ magic_vector + C
end

function get_oxygen_coords(
    backbone::Backbone{3, T},
) where T <: Real
    L = length(backbone)

    CAs = alphacarbons(backbone)
    Cs = carbons(backbone)
    next_Ns = view(nitrogens(backbone), 2:L)

    oxygen_coords = zeros(T, 3, L-1)
    for (i, (CA, C, next_N)) in enumerate(zip(CAs, Cs, next_Ns))
        oxygen_coords[:, i] = get_oxygen_coords(CA, C, next_N)
    end
    
    return oxygen_coords
end

function backbone_with_oxygen(backbone::Backbone{3, T}) where T <: Real
    L = length(backbone)
    oxygen_coords = get_oxygen_coords(backbone)
    coords = zeros(T, 3, 4, L-1)
    coords[:, 1:3, :] = backbone.coords[:, :, 1:L-1]
    coords[:, 4, :] = oxygen_coords
    return Backbone(coords)
end