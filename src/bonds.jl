import Rotations

export get_atom_displacements
export get_atom_distances

export get_bond_vectors
export get_bond_lengths
export get_bond_angles
export get_dihedrals

export ChainedBonds

export append_bonds!
export append_bonds

column_norms(vectors::AbstractMatrix) = reshape(mapslices(norm, vectors, dims=1), :)

function get_atom_displacements(
    backbone::Backbone, start::Integer, step::Integer, stride::Integer,
)
    a1 = backbone[start:stride:end-step].coords
    a2 = backbone[start+step:stride:end].coords
    @assert size(a1, 2) == size(a2, 2)
    displacements = a2 .- a1
    return displacements
end

function get_atom_distances(
    backbone::Backbone, start::Integer, step::Integer, stride::Integer,
)
    displacements = get_atom_displacements(backbone, start, step, stride)
    distances = column_norms(displacements)
    return distances
end


get_bond_vectors(backbone::Backbone) = get_atom_displacements(backbone, 1, 1, 1)

get_bond_lengths(bond_vectors::AbstractMatrix{<:Real}) = column_norms(bond_vectors)
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

# Example
```jldoctest
julia> protein = readpdb("test/data/1ZAK.pdb")
2-element Vector{Chain}:
 Chain A with 220 residues
 Chain B with 220 residues

julia> chain = protein["A"]
Chain A with 220 residues

julia> oxygen_coords(chain) # returns the estimated position of oxygen atoms (~0.05 Å mean deviation)
3×220 Matrix{Float32}:
 22.6697  25.1719  24.7761  25.8559  …  24.7911   22.7649   22.6578   21.24
 15.7257  13.505   13.5151  11.478      15.0888   12.2361   15.8825   14.2933
 21.4295  19.5663  22.8638  25.3283      7.95346   4.81901   3.24164  -0.742424        
```

!!! note
    The `oxygen_coords` function adds oxygen atoms to the backbone using idealized geometry, and oxygens atom will on average deviate [0.05 Å](https://github.com/MurrellGroup/Backboner.jl/blob/main/test/protein/oxygen.jl) from the original positions.
    Moreover, the last oxygen atom is essentially given a random (although deterministic) orientation, as that information is lost when the backbone is reduced to 3 atoms, and there's no next nitrogen atom to compare with.
"""
struct ChainedBonds{T <: Real}
    lengths::Vector{T}
    angles::Vector{T}
    dihedrals::Vector{T}

    function ChainedBonds(lengths::Vector{T}, angles::Vector{T}, dihedrals::Vector{T}) where T
        @assert length(lengths) == length(angles) + 1 == length(dihedrals) + 2
        return new{T}(lengths, angles, dihedrals)
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

    ChainedBonds(backbone::Backbone) = ChainedBonds(get_bond_vectors(backbone))
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
            bond_angle_rotation = Rotations.AngleAxis(π - bonds.angles[1], N...)
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
        bond_angle_rotation = Rotations.AngleAxis(π - bonds.angles[i-2], N...)
        CD_rot1 = bond_angle_rotation * CD_init

        dihedral_rotation = Rotations.AngleAxis(bonds.dihedrals[i-3], n_BC...)
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


get_skip_length(L1, L2, θ) = sqrt(L1^2 + L2^2 - 2*L1*L2*cos(θ))

function get_skip_lengths(lengths::V, angles::V) where V <: AbstractVector{<:Real}
    @assert length(lengths) == length(angles) + 1
    L1s = lengths
    L2s = Iterators.drop(lengths, 1)
    skip_lengths = similar(angles)
    for (L1, L2, (i, θ)) in zip(L1s, L2s, enumerate(angles))
        skip_lengths[i] = get_skip_length(L1, L2, θ)
    end
    return skip_lengths
end

get_skip_lengths(bonds::ChainedBonds{T}) where T = get_skip_lengths(bonds.lengths, bonds.angles)
