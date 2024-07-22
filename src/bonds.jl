using LinearAlgebra
using Rotations: AngleAxis

norms(A::AbstractArray{<:Real}; dims=1) = sqrt.(sum(abs2, A; dims))
dots(A1::AbstractArray{T}, A2::AbstractMatrix{T}; dims=1) where T <: Real = sum(A1 .* A2; dims)
normalize_slices(A::AbstractArray{<:Real}; dims=1) = A ./ norms(A; dims)

get_atom_displacements(backbone::Backbone, start::Integer, step::Integer, stride::Integer) =
    backbone.coords[:, start+step:stride:end] .- backbone.coords[:, start:stride:end-step]

get_atom_distances(backbone::Backbone, start::Integer, step::Integer, stride::Integer) =
    norms(get_atom_displacements(backbone, start, step, stride))

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

function _get_dihedrals(bond_vectors::AbstractMatrix{<:Real})
    crosses = batched_cross(bond_vectors[:, begin:end-1], bond_vectors[:, begin+1:end])
    normalized_crosses = normalize_slices(crosses)
    cross_crosses = batched_cross(normalized_crosses[:, begin:end-1], normalized_crosses[:, begin+1:end])
    normalized_bond_vectors = normalize_slices(bond_vectors)
    sin_values = dots(cross_crosses, normalized_bond_vectors[:, begin+1:end-1])
    cos_values = dots(normalized_crosses[:, begin:end-1], normalized_crosses[:, begin+1:end])
    dihedrals = atan.(sin_values, cos_values)
    return dihedrals
end

get_bond_vectors(backbone::Backbone) = get_atom_displacements(backbone, 1, 1, 1)
get_bond_lengths(backbone::Backbone) = _get_bond_lengths(get_bond_vectors(backbone)) |> vec
get_bond_angles(backbone::Backbone) = _get_bond_angles(get_bond_vectors(backbone)) |> vec
get_dihedrals(backbone::Backbone) = _get_dihedrals(get_bond_vectors(backbone)) |> vec

"""
    ChainedBonds{T <: Real, V <: AbstractVector{T}}

A lazy way to store a backbone as a series of bond lengths, angles, and dihedrals.

# Examples

```jldoctest
julia> backbone = Protein.readpdb("test/data/1ZAK.pdb")["A"].backbone
660-element Backbone{Float64, Matrix{Float64}}:
 [22.346, 17.547, 23.294]
 [22.901, 18.031, 21.993]
 [23.227, 16.793, 21.163]
 [24.115, 16.923, 20.175]
 [24.478, 15.779, 19.336]
 ⋮
 [21.480, 14.668, 4.974]
 [22.041, 14.866, 3.569]
 [21.808, 13.861, 2.734]
 [22.263, 13.862, 1.355]
 [21.085, 14.233, 0.446]

julia> bonds = ChainedBonds(backbone)
ChainedBonds{Float64, Vector{Float64}} with 659 bonds, 658 angles, and 657 dihedrals
```
"""
struct ChainedBonds{T <: Real, V <: AbstractVector{T}}
    lengths::V
    angles::V
    dihedrals::V

    function ChainedBonds{T, V}(lengths::V, angles::V, dihedrals::V) where {T <: Real, V <: AbstractVector{T}}
        length(lengths) == length(angles) + 1 == length(dihedrals) + 2 || 1 >= length(lengths) >= length(angles) == length(dihedrals) == 0 || throw(ArgumentError("lengths, angles, and dihedrals must have compatible lengths"))
        return new{T, V}(lengths, angles, dihedrals)
    end
end

@inline ChainedBonds{T}(lengths::V, angles::V, dihedrals::V) where {T <: Real, V <: AbstractVector{T}} = ChainedBonds{T, V}(lengths, angles, dihedrals)
@inline ChainedBonds(lengths::V, angles::V, dihedrals::V) where {T <: Real, V <: AbstractVector{T}} = ChainedBonds{T, V}(lengths, angles, dihedrals)

get_bond_lengths(bonds::ChainedBonds) = bonds.lengths
get_bond_angles(bonds::ChainedBonds) = bonds.angles
get_bond_dihedrals(bonds::ChainedBonds) = bonds.dihedrals

function Base.reverse!(bonds::ChainedBonds)
    reverse!(bonds.lengths)
    reverse!(bonds.angles)
    reverse!(bonds.dihedrals)
    return bonds
end

Base.reverse(bonds::ChainedBonds) = reverse!(deepcopy(bonds))

function ChainedBonds(backbone::Backbone)
    bond_vectors = get_bond_vectors(backbone)
    lengths = _get_bond_lengths(bond_vectors) |> vec
    angles = _get_bond_angles(bond_vectors) |> vec
    dihedrals = _get_dihedrals(bond_vectors) |> vec
    return ChainedBonds(lengths, angles, dihedrals)
end

Base.length(bonds::ChainedBonds) = length(bonds.lengths)

function Base.show(io::IO, bonds::ChainedBonds)
    print(io, "$(summary(bonds)) with $(length(bonds.lengths)) bonds, $(length(bonds.angles)) angles, and $(length(bonds.dihedrals)) dihedrals")
end


function get_backbone_start(bonds::ChainedBonds{T}) where T
    L = length(bonds) + 1
    l = min(3, L)
    coords = Matrix{T}(undef, 3, l)
    coords[:, 1] = [0, 0, 0]
    if l >= 2
        coords[:, 2] = [bonds.lengths[1], 0, 0]
        if l == 3
            N = normalize([0, 0, 1])
            bond_angle_rotation = AngleAxis(π - bonds.angles[1], N...)
            coords[:, 3] = coords[:, 2] + bond_angle_rotation * [bonds.lengths[2], 0, 0]
        end
    end
    return Backbone(coords)
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
        A, B, C = eachcol(coords[:, i-3:i-1])
        AB, BC = B - A, C - B
        n_AB, n_BC = normalize(AB), normalize(BC)
        CD_init = n_BC * bonds.lengths[i-1]

        N = normalize(cross(n_AB, n_BC)) # wont work if AB == BC        
        bond_angle_rotation = AngleAxis(π - bonds.angles[i-2], N...)
        CD_rot1 = bond_angle_rotation * CD_init

        dihedral_rotation = AngleAxis(bonds.dihedrals[i-3], n_BC...)
        CD_rot2 = dihedral_rotation * CD_rot1

        D = C + CD_rot2
        coords[:, i] = D
    end
    return Backbone(coords)
end


function append_bonds!(bonds::ChainedBonds, lengths::AbstractVector{<:Real}, angles::AbstractVector{<:Real}, dihedrals::AbstractVector{<:Real})
    length(bonds) >= 2 || throw(ArgumentError("bonds must have at least 2 bonds"))
    length(lengths) == length(angles) == length(dihedrals) || throw(ArgumentError("lengths, angles, and dihedrals must have the same length"))
    append!(bonds.lengths, lengths)
    append!(bonds.angles, angles)
    append!(bonds.dihedrals, dihedrals)
    return bonds
end

function append_bonds(bonds::ChainedBonds, lengths::AbstractVector{<:Real}, angles::AbstractVector{<:Real}, dihedrals::AbstractVector{<:Real})
    deepcopy(bonds)
    append_bonds!(bonds, lengths, angles, dihedrals)
    return bonds
end

"""
    append_bonds(backbone, lengths, angles, dihedrals)
"""
function append_bonds(backbone::Backbone, lengths::AbstractVector{<:Real}, angles::AbstractVector{<:Real}, dihedrals::AbstractVector{<:Real})
    length(backbone) >= 3 || throw(ArgumentError("backbone must have at least 3 atoms"))
    length(lengths) == length(angles) == length(dihedrals) || throw(ArgumentError("lengths, angles, and dihedrals must have the same length"))
    last_three_atoms_backbone = backbone[end-2:end]
    bonds_end = ChainedBonds(last_three_atoms_backbone)
    append_bonds!(bonds_end, lengths, angles, dihedrals)
    backbone_end = Backbone(bonds_end, backbone_start=last_three_atoms_backbone)
    new_backbone = Backbone(cat(backbone.coords, backbone_end[4:end].coords, dims=2))
    return new_backbone
end

function prepend_bonds!(bonds::ChainedBonds, lengths::AbstractVector{<:Real}, angles::AbstractVector{<:Real}, dihedrals::AbstractVector{<:Real})
    length(bonds) >= 2 || throw(ArgumentError("bonds must have at least 2 bonds"))
    length(lengths) == length(angles) == length(dihedrals) || throw(ArgumentError("lengths, angles, and dihedrals must have the same length"))
    prepend!(bonds.lengths, lengths)
    prepend!(bonds.angles, angles)
    prepend!(bonds.dihedrals, dihedrals)
    return bonds
end

function prepend_bonds(bonds::ChainedBonds, lengths::AbstractVector{<:Real}, angles::AbstractVector{<:Real}, dihedrals::AbstractVector{<:Real})
    deepcopy(bonds)
    prepend_bonds!(bonds, lengths, angles, dihedrals)
    return bonds
end

"""
    prepend_bonds(backbone, lengths, angles, dihedrals)
"""
function prepend_bonds(backbone::Backbone, lengths::AbstractVector{<:Real}, angles::AbstractVector{<:Real}, dihedrals::AbstractVector{<:Real})
    length(backbone) >= 3 || throw(ArgumentError("backbone must have at least 3 atoms"))
    first_three_atoms_backbone = backbone[1:3]
    bonds_start = ChainedBonds(first_three_atoms_backbone)
    prepend_bonds!(bonds_start, lengths, angles, dihedrals)
    reverse!(bonds_start)
    backbone_start = Backbone(bonds_start, backbone_start=first_three_atoms_backbone[3:-1:1])
    new_backbone = Backbone(cat(backbone_start[end:-1:4].coords, backbone.coords, dims=2))
    return new_backbone
end