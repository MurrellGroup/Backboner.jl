using LinearAlgebra
import Rotations: AngleAxis
import Distances: euclidean

export get_atom_displacements
export get_atom_distances

export get_bond_vectors
export get_bond_lengths
export get_bond_angles
export get_dihedrals

export ChainedBonds

export append_bonds!
export append_bonds

_column_norms(columns::AbstractMatrix) = reshape(mapslices(norm, columns, dims=1), :)

function _pairwise_column_displacements(
    columns1::M, columns2::M
) where M <: AbstractMatrix{<:Real}
    @assert size(columns1) == size(columns2)
    displacements = columns2 .- columns1
    return displacements
end

function _pairwise_column_distances(
    columns1::M, columns2::M, f = euclidean
) where {T <: Real, M <: AbstractMatrix{T}}
    @assert size(columns1, 2) == size(columns2, 2)
    distances = Vector{T}(undef, size(columns1, 2))
    @inbounds for (i, (col1, col2)) in enumerate(zip(eachcol(columns1), eachcol(columns2)))
        distances[i] = f(col1, col2)
    end
    return distances
end

function get_atom_displacements(
    backbone::Backbone, start::Integer, step::Integer, stride::Integer,
)
    displacements = _pairwise_column_displacements(
        @view(backbone.coords[:, start:stride:end-step]),
        @view(backbone.coords[:, start+step:stride:end]))
    return displacements
end

function get_atom_distances(
    backbone::Backbone, start::Integer, step::Integer, stride::Integer,
)
    distances = _pairwise_column_distances(
        @view(backbone.coords[:, start:stride:end-step]),
        @view(backbone.coords[:, start+step:stride:end]))
    return distances
end


get_bond_vectors(backbone::Backbone) = get_atom_displacements(backbone, 1, 1, 1)

get_bond_lengths(bond_vectors::AbstractMatrix{<:Real}) = _column_norms(bond_vectors)
get_bond_lengths(backbone::Backbone) = get_bond_lengths(get_bond_vectors(backbone))

function get_bond_angle(v1::V, v2::V) where V <: AbstractVector{<:Real}
    return acos((-dot(v1, v2)) / (norm(v1) * norm(v2)))
end

function get_bond_angles(bond_vectors::AbstractMatrix{<:Real})
    bond_vector_pairs = zip(eachcol(bond_vectors), Iterators.drop(eachcol(bond_vectors), 1))
    bond_angles = [get_bond_angle(v1, v2) for (v1, v2) in bond_vector_pairs]
    return bond_angles
end

get_bond_angles(backbone::Backbone) = get_bond_angles(get_bond_vectors(backbone))

# source: https://en.wikipedia.org/w/index.php?title=Dihedral_angle&oldid=1182848974#In_polymer_physics
function dihedral_angle(u1::V, u2::V, u3::V) where V <: AbstractVector{<:Real}
    c12, c23 = cross(u1, u2), cross(u2, u3)
    return atan(dot(u2, cross(c12, c23)), norm(u2) * dot(c12, c23))
end

function get_dihedrals(vectors::AbstractMatrix{T}) where T
    dihedrals = Vector{T}(undef, size(vectors, 2) - 2)
    u1s = eachcol(vectors)
    u2s = Iterators.drop(u1s, 1)
    u3s = Iterators.drop(u1s, 2)
    for (i, u1, u2, u3) in zip(eachindex(dihedrals), u1s, u2s, u3s)
        dihedrals[i] = dihedral_angle(u1, u2, u3)
    end
    return dihedrals
end

get_dihedrals(backbone::Backbone) = get_dihedrals(get_bond_vectors(backbone))


"""
    ChainedBonds{T <: Real}

A lazy way to store a backbone as a series of bond lengths, angles, and dihedrals.
It can be instantiated from a Backbone or a matrix of bond vectors.
It can also be used to instantiate a Backbone using the `Backbone(bonds::ChainedBonds)` constructor.

# Examples

```jldoctest
julia> backbone = readpdb("test/data/1ZAK.pdb")["A"].backbone
660-element Backbone{Float32}:
 [22.346, 17.547, 23.294]
 [22.901, 18.031, 21.993]
 [23.227, 16.793, 21.163]
 [24.115, 16.923, 20.175]
 ⋮
 [22.041, 14.866, 3.569]
 [21.808, 13.861, 2.734]
 [22.263, 13.862, 1.355]
 [21.085, 14.233, 0.446]

julia> bonds = ChainedBonds(backbone)
ChainedBonds{Float32} with 659 bonds, 658 angles, and 657 dihedrals
```
"""
struct ChainedBonds{T <: Real}
    lengths::Vector{T}
    angles::Vector{T}
    dihedrals::Vector{T}

    function ChainedBonds(lengths::Vector{T}, angles::Vector{T}, dihedrals::Vector{T}) where T
        @assert length(lengths) == length(angles) + 1 == length(dihedrals) + 2 || 1 >= length(lengths) >= length(angles) == length(dihedrals) == 0
        return new{T}(lengths, angles, dihedrals)
    end
end

function ChainedBonds(lengths::AbstractVector{T}, angles::AbstractVector{T}, dihedrals::AbstractVector{T}) where T
    return ChainedBonds(Vector(lengths), Vector(angles), Vector(dihedrals))
end

function ChainedBonds(vectors::AbstractMatrix{<:Real})
    lengths = get_bond_lengths(vectors)
    angles = get_bond_angles(vectors)
    dihedrals = get_dihedrals(vectors)
    return ChainedBonds(lengths, angles, dihedrals)
end

function ChainedBonds(backbone::Backbone)
    return ChainedBonds(get_bond_vectors(backbone))
end

Base.:(==)(b1::ChainedBonds, b2::ChainedBonds) = b1.lengths == b2.lengths && b1.angles == b2.angles && b1.dihedrals == b2.dihedrals
Base.:(≈)(b1::ChainedBonds, b2::ChainedBonds) = b1.lengths ≈ b2.lengths && b1.angles ≈ b2.angles && b1.dihedrals ≈ b2.dihedrals
Base.length(bonds::ChainedBonds) = length(bonds.lengths)
Base.size(bonds::ChainedBonds) = Tuple(length(bonds))

Base.summary(bonds::ChainedBonds) = string(typeof(bonds))

function Base.show(io::IO, bonds::ChainedBonds)
    print(io, "$(summary(bonds)) with $(length(bonds.lengths)) bonds, $(length(bonds.angles)) angles, and $(length(bonds.dihedrals)) dihedrals")
end


function get_first_points(bonds::ChainedBonds{T}) where T
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

    return coords
end

# first_points don't get adjusted to fit the bonds
# first_points needs to have at least 3 columns
function Backbone(
    bonds::ChainedBonds{T};
    first_points::AbstractMatrix{<:Real} = get_first_points(bonds),
) where T
    L = length(bonds) + 1
    coords = fill(T(NaN), 3, L)

    l = size(first_points, 2)
    @assert 3 <= l <= L
    coords[:, 1:l] = first_points

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

    backbone = Backbone(coords)
    return backbone
end


function append_bonds!(
    bonds::ChainedBonds,
    lengths::AbstractVector{<:Real},
    angles::AbstractVector{<:Real},
    dihedrals::AbstractVector{<:Real},
)
    @assert length(bonds) >= 2
    @assert length(lengths) == length(angles) == length(dihedrals)
    append!(bonds.lengths, lengths)
    append!(bonds.angles, angles)
    append!(bonds.dihedrals, dihedrals)
    return bonds
end

function append_bonds(
    bonds::ChainedBonds,
    lengths::AbstractVector{<:Real},
    angles::AbstractVector{<:Real},
    dihedrals::AbstractVector{<:Real},
)
    deepcopy(bonds)
    append_bonds!(bonds, lengths, angles, dihedrals)
    return bonds
end

"""
    append_bonds(backbone, lengths, angles, dihedrals)
"""
function append_bonds(
    backbone::Backbone,
    lengths::AbstractVector{<:Real},
    angles::AbstractVector{<:Real},
    dihedrals::AbstractVector{<:Real},
)
    @assert length(backbone) >= 3
    @assert length(lengths) == length(angles) == length(dihedrals)
    last_three_atoms_backbone = backbone[end-2:end]
    bonds_end = ChainedBonds(last_three_atoms_backbone)
    append_bonds!(bonds_end, lengths, angles, dihedrals)
    backbone_end = Backbone(bonds_end, first_points = last_three_atoms_backbone)
    new_backbone = Backbone(cat(backbone.coords, backbone_end[4:end].coords, dims=2))
    return new_backbone
end