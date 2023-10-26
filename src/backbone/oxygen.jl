const magic_vector = [-0.672, -1.034, 0.003]

function get_oxygen_coords(
    CA::AbstractVector{T}, C::AbstractVector{T}, next_N::AbstractVector{T},
) where T <: Real
    R = get_rotation_matrix(CA, C, next_N)
    return R' \ magic_vector + C
end

function get_oxygen_coords(
    backbone::Backbone{T},
) where T <: Real
    L = nresidues(backbone)

    CAs = alphacarbons(backbone)
    Cs = carbons(backbone)
    next_Ns = view(nitrogens(backbone), 2:L)

    oxygen_coords = zeros(T, L-1, 3)
    for (i, (CA, C, next_N)) in enumerate(zip(CAs, Cs, next_Ns))
        oxygen_coords[i, :] = get_oxygen_coords(CA, C, next_N)
    end
    
    return oxygen_coords
end

function ResidueArray{4}(backbone::Backbone{T}) where T <: Real
    L = nresidues(backbone)
    oxygen_coords = get_oxygen_coords(backbone)
    coords = zeros(T, L-1, 4, 3)
    coords[:, 1:3, :] = backbone.coords[1:L-1, :, :]
    coords[:, 4, :] = oxygen_coords
    return BackboneAndOxygen(coords)
end

oxygen_coord_matrix(ra::ResidueArray) = _ith_column_positions(ra, 4)
oxygens(ra::ResidueArray) = eachrow(_ith_column_positions(ra, 4))