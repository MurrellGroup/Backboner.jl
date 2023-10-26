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