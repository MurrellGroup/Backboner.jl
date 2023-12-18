export get_atom_displacements, get_atom_distances
export get_bond_vectors, get_bond_lengths
export get_bond_angles, get_dihedrals
export ChainedBonds

column_norms(vectors::AbstractMatrix) = reshape(mapslices(norm, vectors, dims=1), :)

function get_atom_displacements(
    backbone::Backbone{A}, atom1::Integer, atom2::Integer, residue_offset::Integer = 0,
) where A
    @assert 1 <= atom1 <= A && 1 <= atom2 <= A "Backbone{$N} does not have atoms $atom1 and $atom2"
    coords = backbone.coords
    displacements = @view(coords[:, atom2, 1+residue_offset:end]) .- @view(coords[:, atom1, 1:end-residue_offset])
    return displacements
end

function get_atom_distances(
    backbone::Backbone{A}, atom1::Integer, atom2::Integer, residue_offset::Integer = 0,
) where A
    @assert 1 <= atom1 <= A && 1 <= atom2 <= A "Backbone{$N} does not have atoms $atom1 and $atom2"
    displacements = get_atom_displacements(backbone, atom1, atom2, residue_offset)
    distances = column_norms(displacements)
    return distances
end


function get_bond_vectors(backbone::Backbone{A}) where A
    backbone1 = Backbone{1}(backbone)
    bond_vectors = get_atom_displacements(backbone1, 1, 1, 1)
    return bond_vectors
end


function get_bond_lengths(bond_vectors::AbstractMatrix{T}) where T
    bond_lengths = column_norms(bond_vectors)
    return bond_lengths
end

get_bond_lengths(backbone::Backbone) = get_bond_lengths(get_bond_vectors(backbone))


get_bond_angle(v1::V, v2::V) where V <: AbstractVector{<:Real} = acos((-dot(v1, v2)) / (norm(v1) * norm(v2)))

function get_bond_angles(bond_vectors::AbstractMatrix{T}) where T
    bond_vector_pairs = zip(eachcol(bond_vectors), Iterators.drop(eachcol(bond_vectors), 1))
    bond_angles = [get_bond_angle(v1, v2) for (v1, v2) in bond_vector_pairs]
    return bond_angles
end

get_bond_angles(backbone::Backbone) = get_bond_angles(get_bond_vectors(backbone))

# source: en.wikipedia.org/wiki/Dihedral_angle#In_polymer_physics
function dihedral_angle(u1::V, u2::V, u3::V) where V <: AbstractVector{<:Real}
    c12, c23 = cross(u1, u2), cross(u2, u3)
    atan(dot(u2, cross(c12, c23)), norm(u2) * dot(c12, c23))
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

Can be instantiated from a Backbone or a matrix of bond vectors.

Can be used to instantiate a Backbone using the `Backbone{A}(bonds::ChainedBonds)` constructor.
"""
struct ChainedBonds{T <: Real}
    lengths::AbstractVector{T}
    angles::AbstractVector{T}
    dihedrals::AbstractVector{T}

    function ChainedBonds(vectors::AbstractMatrix{T}) where T
        lengths = get_bond_lengths(vectors)
        angles = get_bond_angles(vectors)
        dihedrals = get_dihedrals(vectors)

        return new{T}(lengths, angles, dihedrals)
    end

    function ChainedBonds(backbone::Backbone{A, T}) where {A, T}
        return ChainedBonds(get_bond_vectors(backbone))
    end
end

Base.:(==)(b1::ChainedBonds, b2::ChainedBonds) = b1.lengths == b2.lengths && b1.angles == b2.angles && b1.dihedrals == b2.dihedrals
Base.:(≈)(b1::ChainedBonds, b2::ChainedBonds) = b1.lengths ≈ b2.lengths && b1.angles ≈ b2.angles && b1.dihedrals ≈ b2.dihedrals
Base.length(bonds::ChainedBonds) = length(bonds.lengths)
Base.size(bonds::ChainedBonds) = Tuple(length(bonds))

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

# first points currently don't get adjusted to fit the bonds
function Backbone{ATOMS_PER_RESIDUE}(
    bonds::ChainedBonds{T};
    first_points::AbstractMatrix{T} = get_first_points(bonds),
) where {ATOMS_PER_RESIDUE, T}
    @assert (length(bonds) + 1) % ATOMS_PER_RESIDUE == 0 "Invalid number of atoms per residue in backbone"

    L = length(bonds) + 1
    coords = fill(T(NaN), 3, L)

    l = size(first_points, 2)
    @assert l <= L
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

    backbone = Backbone(reshape(coords, 3, ATOMS_PER_RESIDUE, :))
    return backbone
end

Backbone(bonds::ChainedBonds; kwargs...) = Backbone{1}(bonds; kwargs...)

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
