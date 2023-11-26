export locs_and_rots_to_backbone, backbone_to_locs_and_rots

const STANDARD_TRIANGLE_ANGSTROM = Float32[
    -0.876  -0.248   1.124;
     0.656  -0.656   0.000;
     0.000   0.000   0.000;
] #  N       Ca      C

const STANDARD_TRIANGLE_NM = STANDARD_TRIANGLE_ANGSTROM .* 0.1

@inline standardtriangle(location::AbstractVector, quatrots::Rotations.QuatRotation) = quatrots * STANDARD_TRIANGLE_ANGSTROM .+ location
@inline standardtriangle(location::AbstractVector, rot_matrix::AbstractMatrix) = rot_matrix * STANDARD_TRIANGLE_ANGSTROM .+ location

# 3x3xL array of rotation matrices
"""
    locs_and_rots_to_backbone(locations, rot_matrices; unit=:angstrom)

Returns a backbone with the given locations and rotation matrices of residues.
If unit is :nm, the locations are converted to angstroms by multiplying them by 10.
"""
function locs_and_rots_to_backbone(locations::AbstractMatrix{T}, rot_matrices::AbstractArray{T, 3}; unit::Symbol=:angstrom) where T
    unit == :nm && (locations *= 10.0)
    @assert size(locations, 1) == 3 "locations must be of size 3xL"
    @assert size(rot_matrices, 1) == size(rot_matrices, 2) == 3 "rot_matrices must be a 3x3xL array"
    @assert size(locations, 2) == size(rot_matrices, 3) "The second dimension of locations must be the same as the third dimension of rot_matrices"
    triangles = standardtriangle.(eachcol(locations), eachslice(rot_matrices, dims=3))
    return Backbone(stack(triangles))
end

# length L vector of Rotations.QuatRotation
function locs_and_rots_to_backbone(locations::AbstractMatrix{T}, quatrots::AbstractVector{Rotations.QuatRotation{T}}; unit::Symbol=:angstrom) where T
    unit == :nm && (locations *= 10.0)
    @assert size(locations, 1) == 3 "locations must be of size 3xL"
    @assert size(locations, 2) == length(quatrots) "The second dimension of locations must be the same as the length of quatrots"
    triangles = standardtriangle.(eachcol(locations), quatrots)
    return Backbone(stack(triangles))
end

# 4xL matrix of quaternions
function locs_and_rots_to_backbone(locations::AbstractMatrix{T}, quatrot_matrix::AbstractMatrix{T}; unit::Symbol=:angstrom) where T
    @assert size(quatrot_matrix, 1) == 4 "quatrot_matrix must be a 4xL array"
    quatrots = Rotations.QuatRotation.(eachcol(quatrot_matrix))
    return locs_and_rots_to_backbone(locations, quatrots, unit=unit)
end

"""
    backbone_to_locs_and_rots(backbone, unit=:angstrom)

Returns the locations and rotation matrices of residues in a backbone,
according to a defined standard triangle (`Backboner.STANDARD_TRIANGLE_ANGSTROM`).
"""
function backbone_to_locs_and_rots(backbone::Backbone{N, T}, unit::Symbol=:angstrom) where {N, T}
    @assert 3 <= N <= 4 "backbone must have 3 or 4 atoms per residue"
    L = length(backbone)
    locations = Array{T}(undef, 3, L)
    rot_matrices = Array{T}(undef, 3, 3, L)
    for i in 1:L
        triangle = view(backbone, :, 1:3, i)
        location = reduce(+, triangle, dims=2) / 3
        rot_matrix = (triangle .- location) * pinv(STANDARD_TRIANGLE_ANGSTROM)
        locations[:, i] = location
        rot_matrices[:, :, i] = rot_matrix
    end
    unit == :nm && (locations *= 0.1)
    return locations, rot_matrices
end