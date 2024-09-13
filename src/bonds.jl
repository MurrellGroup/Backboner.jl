# TODO: see if views would be differentiable (GPU?)

"""
    get_displacements(backbone, start, step, stride)

Get the displacements between points in `backbone`, where `start` is the index of the first point,
`step` is the offset of the other points along the backbone (e.g. `1` for bond vectors),
and `stride` is the number of points to skip after each step.
"""
get_displacements(backbone::Backbone, start::Integer, step::Integer, stride::Integer) =
    backbone.coords[:, start+step:stride:end] .- backbone.coords[:, start:stride:end-step]

"""
    get_distances(backbone, start, step, stride)

Get the Euclidean distances between points in `backbone`, where `start` is the index of the first point,
`step` is the offset of the other points along the backbone (e.g. `1` for bond lengths),
and `stride` is the number of points to skip after each step.
"""
get_distances(backbone::Backbone, start::Integer, step::Integer, stride::Integer) =
    norms(get_displacements(backbone, start, step, stride))

_get_bond_lengths(bond_vectors::AbstractMatrix{<:Real}) = norms(bond_vectors)

function _get_bond_angles(bond_vectors::AbstractMatrix{T}) where T
    us = bond_vectors[:, begin:end-1]
    vs = bond_vectors[:, begin+1:end]
    return π .- acos.(clamp.(dots(us, vs) ./ (norms(us) .* norms(vs)), -one(T), one(T)))
end

function batched_cross(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T <: Real
    C1 = selectdim(A, 1, 2) .* selectdim(B, 1, 3) .- selectdim(A, 1, 3) .* selectdim(B, 1, 2)
    C2 = selectdim(A, 1, 3) .* selectdim(B, 1, 1) .- selectdim(A, 1, 1) .* selectdim(B, 1, 3)
    C3 = selectdim(A, 1, 1) .* selectdim(B, 1, 2) .- selectdim(A, 1, 2) .* selectdim(B, 1, 1)
    return [C1 C2 C3]'
end

function _get_torsion_angles(bond_vectors::AbstractMatrix{<:Real})
    crosses = batched_cross(bond_vectors[:, begin:end-1], bond_vectors[:, begin+1:end])
    normalized_crosses = normalize_slices(crosses)
    cross_crosses = batched_cross(normalized_crosses[:, begin:end-1], normalized_crosses[:, begin+1:end])
    normalized_bond_vectors = normalize_slices(bond_vectors)
    sin_values = dots(cross_crosses, normalized_bond_vectors[:, begin+1:end-1])
    cos_values = dots(normalized_crosses[:, begin:end-1], normalized_crosses[:, begin+1:end])
    torsion_angles = atan.(sin_values, cos_values)
    return torsion_angles
end

get_bond_vectors(backbone::Backbone) = get_displacements(backbone, 1, 1, 1)
get_bond_lengths(backbone::Backbone) = _get_bond_lengths(get_bond_vectors(backbone)) |> vec
get_bond_angles(backbone::Backbone) = _get_bond_angles(get_bond_vectors(backbone)) |> vec
get_torsion_angles(backbone::Backbone) = _get_torsion_angles(get_bond_vectors(backbone)) |> vec

"""
    ChainedBonds{T<:Real,V<:AbstractVector{T}}

A lazy way to store a backbone as a series of bond lengths, bond angles, and torsion_angles.

# Examples

```jldoctest
julia> backbone = Backbone(randn(Float32, 3, 660));

julia> bonds = ChainedBonds(backbone)
ChainedBonds{Float64, Vector{Float64}} with 659 bond lengths, 658 bond angles, and 657 torsion angles
```
"""
struct ChainedBonds{T<:Real,V<:AbstractVector{T}}
    bond_lengths::V
    bond_angles::V
    torsion_angles::V

    function ChainedBonds{T,V}(bond_lengths::V, bond_angles::V, torsion_angles::V) where {T<:Real,V<:AbstractVector{T}}
        length(bond_lengths) == length(bond_angles) + 1 == length(torsion_angles) + 2 || 1 >= length(bond_lengths) >= length(bond_angles) == length(torsion_angles) == 0 ||
            throw(ArgumentError("bond_lengths, bond_angles, and torsion_angles must have compatible lengths"))
        return new{T,V}(bond_lengths, bond_angles, torsion_angles)
    end
end

@inline ChainedBonds{T}(bond_lengths::V, bond_angles::V, torsion_angles::V) where {T<:Real,V<:AbstractVector{T}} = ChainedBonds{T,V}(bond_lengths, bond_angles, torsion_angles)
@inline ChainedBonds(bond_lengths::V, bond_angles::V, torsion_angles::V) where {T<:Real,V<:AbstractVector{T}} = ChainedBonds{T,V}(bond_lengths, bond_angles, torsion_angles)

get_bond_lengths(bonds::ChainedBonds) = bonds.bond_lengths
get_bond_angles(bonds::ChainedBonds) = bonds.bond_angles
get_torsion_angles(bonds::ChainedBonds) = bonds.torsion_angles

function Base.reverse!(bonds::ChainedBonds)
    reverse!(bonds.bond_lengths)
    reverse!(bonds.bond_angles)
    reverse!(bonds.torsion_angles)
    return bonds
end

Base.reverse(bonds::ChainedBonds) = reverse!(deepcopy(bonds))

function ChainedBonds(backbone::Backbone)
    bond_vectors = get_bond_vectors(backbone)
    bond_lengths = _get_bond_lengths(bond_vectors) |> vec
    bond_angles = _get_bond_angles(bond_vectors) |> vec
    torsion_angles = _get_torsion_angles(bond_vectors) |> vec
    return ChainedBonds(bond_lengths, bond_angles, torsion_angles)
end

Base.length(bonds::ChainedBonds) = length(bonds.bond_lengths)

function Base.show(io::IO, bonds::ChainedBonds)
    print(io, "$(summary(bonds)) with $(length(bonds.bond_lengths)) bond lengths, $(length(bonds.bond_angles)) bond angles, and $(length(bonds.torsion_angles)) torsion angles")
end


function get_backbone_start(bonds::ChainedBonds{T}) where T
    L = length(bonds) + 1
    l = min(3, L)
    coords = Matrix{T}(undef, 3, l)
    coords[:, 1] = [0, 0, 0]
    if l >= 2
        coords[:, 2] = [bonds.bond_lengths[1], 0, 0]
        if l == 3
            N = normalize([0, 0, 1])
            bond_angle_rotation = AngleAxis(π - bonds.bond_angles[1], N...)
            coords[:, 3] = coords[:, 2] + bond_angle_rotation * [bonds.bond_lengths[2], 0, 0]
        end
    end
    return Backbone(coords)
end

function get_next_point(bonds, i, A, B, C)
    AB, BC = B - A, C - B
    n_AB, n_BC = normalize(AB), normalize(BC)
    CD_init = n_BC * bonds.bond_lengths[i-1]

    N = normalize(cross(n_AB, n_BC)) # wont work if AB == BC        
    bond_angle_rotation = AngleAxis(π - bonds.bond_angles[i-2], N...)
    CD_rot1 = bond_angle_rotation * CD_init

    dihedral_rotation = AngleAxis(bonds.torsion_angles[i-3], n_BC...)
    CD_rot2 = dihedral_rotation * CD_rot1

    D = C + CD_rot2
    return D
end

# backbone_start don't get adjusted to fit the bonds
function Backbone(
    bonds::ChainedBonds{T};
    backbone_start::Backbone{T} = get_backbone_start(bonds),
) where T
    L = length(bonds) + 1
    coords = Matrix{T}(undef, 3, L)
    l = size(backbone_start.coords, 2)
    @assert 3 <= l <= L
    coords[:, 1:l] = backbone_start.coords
    for i in l+1:L
        coords[:, i] = get_next_point(bonds, i, eachcol(coords[:, i-3:i-1])...)
    end
    return Backbone(coords)
end


function append_bonds!(bonds::ChainedBonds, bond_lengths::AbstractVector{<:Real}, bond_angles::AbstractVector{<:Real}, torsion_angles::AbstractVector{<:Real})
    length(bonds) >= 2 || throw(ArgumentError("bonds must have at least 2 bonds"))
    length(bond_lengths) == length(bond_angles) == length(torsion_angles) || throw(ArgumentError("bond_lengths, bond_angles, and torsion_angles must have the same length"))
    append!(bonds.bond_lengths, bond_lengths)
    append!(bonds.bond_angles, bond_angles)
    append!(bonds.torsion_angles, torsion_angles)
    return bonds
end

@inline append_bonds(bonds::ChainedBonds, args...) = append_bonds!(deepcopy(bonds), args...)

"""
    append_bonds(backbone, bond_lengths, bond_angles, torsion_angles)
"""
function append_bonds(backbone::Backbone, bond_lengths::AbstractVector{<:Real}, bond_angles::AbstractVector{<:Real}, torsion_angles::AbstractVector{<:Real})
    length(backbone) >= 3 || throw(ArgumentError("backbone must have at least 3 points"))
    length(bond_lengths) == length(bond_angles) == length(torsion_angles) || throw(ArgumentError("bond_lengths, bond_angles, and torsion_angles must have the same length"))
    last_three_points_backbone = backbone[end-2:end]
    bonds_end = ChainedBonds(last_three_points_backbone)
    append_bonds!(bonds_end, bond_lengths, bond_angles, torsion_angles)
    backbone_end = Backbone(bonds_end, backbone_start=last_three_points_backbone)
    new_backbone = Backbone(cat(backbone.coords, backbone_end[4:end].coords, dims=2))
    return new_backbone
end

function prepend_bonds!(bonds::ChainedBonds, bond_lengths::AbstractVector{<:Real}, bond_angles::AbstractVector{<:Real}, torsion_angles::AbstractVector{<:Real})
    length(bonds) >= 2 || throw(ArgumentError("bonds must have at least 2 bonds"))
    length(bond_lengths) == length(bond_angles) == length(torsion_angles) || throw(ArgumentError("bond_lengths, bond_angles, and torsion_angles must have the same length"))
    prepend!(bonds.bond_lengths, bond_lengths)
    prepend!(bonds.bond_angles, bond_angles)
    prepend!(bonds.torsion_angles, torsion_angles)
    return bonds
end

@inline prepend_bonds(bonds::ChainedBonds, args...) = prepend_bonds!(deepcopy(bonds), args...)

"""
    prepend_bonds(backbone, bond_lengths, bond_angles, torsion_angles)
"""
function prepend_bonds(backbone::Backbone, bond_lengths::AbstractVector{<:Real}, bond_angles::AbstractVector{<:Real}, torsion_angles::AbstractVector{<:Real})
    length(backbone) >= 3 || throw(ArgumentError("backbone must have at least 3 points"))
    first_three_points_backbone = backbone[1:3]
    bonds_start = ChainedBonds(first_three_points_backbone)
    prepend_bonds!(bonds_start, bond_lengths, bond_angles, torsion_angles)
    reverse!(bonds_start)
    backbone_start = Backbone(bonds_start, backbone_start=first_three_points_backbone[3:-1:1])
    new_backbone = Backbone(cat(backbone_start[end:-1:4].coords, backbone.coords, dims=2))
    return new_backbone
end