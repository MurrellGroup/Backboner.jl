export
    get_atom_displacements,
    get_atom_distances,
    get_bond_vectors,
    get_bond_lengths,
    get_bond_angles,
    get_dihedrals,
    ChainedBonds,
    append_bonds!,
    append_bonds

using LinearAlgebra
import Rotations: AngleAxis

column_sums(columns::AbstractMatrix{<:Real}) = vec(sum(columns, dims=1))
column_norms(columns::AbstractMatrix{<:Real}) = sqrt.(column_sums(abs2.(columns)))
column_dots(columns1::M, columns2::M) where M <: AbstractMatrix{<:Real} = column_sums(columns1 .* columns2)

function get_atom_displacements(backbone::Backbone, start::Integer, step::Integer, stride::Integer)
    return @view(backbone.coords[:, start+step:stride:end]) .- @view(backbone.coords[:, start:stride:end-step])
end

function get_atom_distances(backbone::Backbone, start::Integer, step::Integer, stride::Integer)
    return column_norms(get_atom_displacements(backbone, start, step, stride))
end

get_bond_lengths(bond_vectors::AbstractMatrix{<:Real}) = column_norms(bond_vectors)

calculate_bond_angle(u::V, v::V) where {T <: Real, V <: AbstractVector{T}} = π - acos(clamp(dot(u, v) / (norm(u) * norm(v)), -one(T), one(T)))

function get_bond_angles(bond_vectors::AbstractMatrix{T}) where T <: Real
    us = @view(bond_vectors[:, 1:end-1])
    vs = @view(bond_vectors[:, 2:end])
    return π .- acos.(clamp.(column_dots(us, vs) ./ (column_norms(us) .* column_norms(vs)), -one(T), one(T)))
end

# source: https://en.wikipedia.org/w/index.php?title=Dihedral_angle&oldid=1182848974#In_polymer_physics
calculate_dihedral_angle(u1::V, u2::V, u3::V) where V <: AbstractVector{<:Real} = let c12 = cross(u1, u2), c23 = cross(u2, u3)
    atan(dot(u2, cross(c12, c23)), norm(u2) * dot(c12, c23))
end

function get_dihedrals(bond_vectors::AbstractMatrix{<:Real})
    len_prot = size(bond_vectors, 2) + 1
    lengths = sqrt.(bond_vectors[1,:].^2 .+ bond_vectors[2,:].^2 .+ bond_vectors[3,:].^2)

    crosses  = stack(cross.(eachcol(bond_vectors[:,1:len_prot-2]), eachcol(bond_vectors[:,2:len_prot-1])),dims=2)
    lengths_cross = reshape(sqrt.(crosses[1,:].^2 .+ crosses[2,:].^2 .+ crosses[3,:].^2),1,:)
    normalized_crosses = crosses ./ lengths_cross
    cos_theta = dot.(eachcol(normalized_crosses[:,1:len_prot-3]), eachcol(normalized_crosses[:,2:len_prot-2]))

    cross_crosses = stack(cross.(eachcol(normalized_crosses[:,1:len_prot-3]), eachcol(normalized_crosses[:,2:len_prot-2])),dims=2)
    normalized_bond_vectors = (bond_vectors ./ reshape(lengths,1,:))[:,2:len_prot-2]

    sin_theta = dot.(eachcol(cross_crosses), eachcol(normalized_bond_vectors))

    thetas = atan.(sin_theta, cos_theta)
    return thetas
end

get_bond_vectors(backbone::Backbone) = get_atom_displacements(backbone, 1, 1, 1)
get_bond_lengths(backbone::Backbone) = get_atom_distances(backbone, 1, 1, 1)
get_bond_angles(backbone::Backbone) = get_bond_angles(get_bond_vectors(backbone))
get_dihedrals(backbone::Backbone) = get_dihedrals(get_bond_vectors(backbone))


"""
    ChainedBonds{T <: Real, V <: AbstractVector{T}}

A lazy way to store a backbone as a series of bond lengths, angles, and dihedrals.
It can be instantiated from a Backbone or a matrix of bond vectors.
It can also be used to instantiate a Backbone using the `Backbone(bonds::ChainedBonds)` constructor.

# Examples

```jldoctest
julia> backbone = Protein.readpdb("test/data/1ZAK.pdb")["A"].backbone
3×660 Backbone{Float32, Matrix{Float32}}:
 22.346  22.901  23.227  24.115  24.478  25.289  26.091  26.814  …  23.137  22.572  21.48   22.041  21.808  22.263  21.085
 17.547  18.031  16.793  16.923  15.779  14.65   14.958  13.827     13.041  14.235  14.668  14.866  13.861  13.862  14.233
 23.294  21.993  21.163  20.175  19.336  20.009  21.056  21.652      5.676   5.844   4.974   3.569   2.734   1.355   0.446

julia> bonds = ChainedBonds(backbone)
ChainedBonds{Float32, Vector{Float32}} with 659 bonds, 658 angles, and 657 dihedrals
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

function ChainedBonds(backbone::Backbone)
    bond_vectors = get_bond_vectors(backbone)
    lengths = get_bond_lengths(bond_vectors)
    angles = get_bond_angles(bond_vectors)
    dihedrals = get_dihedrals(bond_vectors)
    return ChainedBonds(lengths, angles, dihedrals)
end

Base.length(bonds::ChainedBonds) = length(bonds.lengths)

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