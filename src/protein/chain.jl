export Chain
export has_assigned_ss
export nitrogen_coords, alphacarbon_coords, carbonyl_coords
export nitrogen_alphacarbon_distances, alphacarbon_carbonyl_distances, carbonyl_nitrogen_distances
export phi_angles, psi_angles, omega_angles

"""
    Chain <: AbstractVector{Residue}

A `Chain` represents a chain of a protein, and is a vector of `Residue`s, which are instantiated from indexing the chain.

## Fields
- `id::AbstractString`: A string identifier (usually a single letter).
- `backbone::Backbone`: A backbone with a length divisible by 3, to ensure 3 atoms per residue (N, Ca, C).
- `modelnum::Int`: The model number of the chain.
- `resnums::Vector{Int}`: storing the residue numbers.
- `aavector::Vector{Char}`: storing the amino acid sequence.
- `ssvector::Vector{Char}`: storing the secondary structure.
"""
struct Chain <: AbstractVector{Residue}
    id::AbstractString
    backbone::Backbone
    modelnum::Int
    resnums::Vector{Int}
    aavector::Vector{Char} # TODO: move over to BioStructures.LongAA (from BioSequences)
    ssvector::Vector{Char}

    function Chain(
        id::AbstractString,
        backbone::Backbone;
        modelnum::Int = 1,
        resnums::Vector{Int} = collect(1:length(backbone) ÷ 3),
        aavector::Vector{Char} = fill('G', length(backbone) ÷ 3),
        ssvector::Union{Vector{Char}, Vector{<:Integer}} = fill(' ', length(backbone) ÷ 3),
    )
        L, r = divrem(length(backbone), 3)
        iszero(r) || throw(ArgumentError("backbone must have a length divisible by 3"))
        length(resnums) == L || throw(ArgumentError("length of resnums must be equal to length of backbone divided by 3"))
        length(aavector) == L || throw(ArgumentError("length of aavector must be equal to length of backbone divided by 3"))
        length(ssvector) == L || throw(ArgumentError("length of ssvector must be equal to length of backbone divided by 3"))
        ssvector isa Vector{<:Integer} && (ssvector = get.(('-', 'H', 'E'), ssvector, ' '))

        return new(id, backbone, modelnum, resnums, aavector, ssvector)
    end

    Chain(backbone::Backbone; kwargs...) = Chain("_", backbone; kwargs...) 
end

@inline Base.:(==)(chain1::Chain, chain2::Chain) = chain1.id == chain2.id && chain1.backbone == chain2.backbone && chain1.ssvector == chain2.ssvector
@inline Base.length(chain::Chain) = length(chain.backbone) ÷ 3
@inline Base.size(chain::Chain) = Tuple(length(chain))
@inline Base.getindex(chain::Chain, i::Integer) = Residue(chain.resnums[i], chain.aavector[i], chain.ssvector[i])

Base.summary(chain::Chain) = "Chain $(chain.id) with $(length(chain)) residue$(length(chain) == 1 ? "" : "s")"
Base.show(io::IO, chain::Chain) = print(io, summary(chain))

@inline Base.getindex(protein::AbstractVector{Chain}, id::AbstractString) = protein[findfirst(chain -> chain.id == id, protein)]
@inline Base.getindex(protein::AbstractVector{Chain}, id::Symbol) = protein[String(id)]

has_assigned_ss(ssvector::Vector{Char}) = all(!=(' '), ssvector)
has_assigned_ss(chain::Chain) = has_assigned_ss(chain.ssvector)
has_assigned_ss(protein::AbstractVector{Chain}) = all(has_assigned_ss, protein)

nitrogen_coords(backbone::Backbone) = (@view backbone[1:3:end]).coords
alphacarbon_coords(backbone::Backbone) = (@view backbone[2:3:end]).coords
carbonyl_coords(backbone::Backbone) = (@view backbone[3:3:end]).coords
# oxygen_coords function in src/protein/oxygen.jl
# TODO: hydrogen_coords

"""
    nitrogen_coords(chain::Chain)

Returns the coordinates of all nitrogen atoms in a chain, as a 3xN matrix.
"""
nitrogen_coords(chain::Chain) = nitrogen_coords(chain.backbone)

"""
    alphacarbon_coords(chain::Chain)

Returns the coordinates of all alphacarbon atoms in a chain, as a 3xN matrix.
"""
alphacarbon_coords(chain::Chain) = alphacarbon_coords(chain.backbone)

"""
    carbonyl_coords(chain::Chain)

Returns the coordinates of all carbonyl atoms in a chain, as a 3xN matrix.
"""
carbonyl_coords(chain::Chain) = carbonyl_coords(chain.backbone)

nitrogen_alphacarbon_distances(backbone::Backbone) = get_atom_distances(backbone, 1, 1, 3)
alphacarbon_carbonyl_distances(backbone::Backbone) = get_atom_distances(backbone, 2, 1, 3)
carbonyl_nitrogen_distances(backbone::Backbone) = get_atom_distances(backbone, 3, 1, 3)

"""
    nitrogen_alphacarbon_distances(chain::Chain)

Calculate the distances between all pairs of contiguous nitrogen and alpha-carbon atoms in a chain.
Returns a vector of distances of length `length(chain)`.
"""
nitrogen_alphacarbon_distances(chain::Chain) = nitrogen_alphacarbon_distances(chain.backbone)

"""
    alphacarbon_carbonyl_distances(chain::Chain)

Calculate the distances between all pairs of contiguous alpha-carbon and carbonyl atoms in a chain.
Returns a vector of distances of length `length(chain)`.
"""
alphacarbon_carbonyl_distances(chain::Chain) = alphacarbon_carbonyl_distances(chain.backbone)

"""
    carbonyl_nitrogen_distances(chain::Chain)

Calculate the distances between all pairs of contiguous carbonyl and nitrogen atoms in a chain.
Returns a vector of distances of length `length(chain) - 1`.
"""
carbonyl_nitrogen_distances(chain::Chain) = carbonyl_nitrogen_distances(chain.backbone)

"""
    phi_angles(chain::Chain)
    phi_angles(bonds::ChainedBonds)

Returns the phi (φ) angles of a chain. A `ChainedBonds` object can be passed instead to avoid recalculating the dihedrals.
"""
phi_angles(chain::Chain) = phi_angles(ChainedBonds(chain))
phi_angles(bonds::ChainedBonds) = @view bonds.dihedrals[3:3:end]

"""
    psi_angles(chain::Chain)
    psi_angles(bonds::ChainedBonds)

Returns the psi (ψ) angles of a chain.
"""
psi_angles(chain::Chain) = psi_angles(ChainedBonds(chain))
psi_angles(bonds::ChainedBonds) = @view bonds.dihedrals[1:3:end]

"""
    omega_angles(chain::Chain)
    omega_angles(bonds::ChainedBonds)

Returns the phi (Ω) angles of a chain.
"""
omega_angles(chain::Chain) = omega_angles(ChainedBonds(chain))
omega_angles(bonds::ChainedBonds) = @view bonds.dihedrals[2:3:end]


# TODO: functions for getting beta-carbon and hydrogen atom positions? is that even possible without knowing AAs and molecular dynamics?

# TODO: append_residues! function. also get the Residue type sorted out. user shouldn't need to create it manually.
#=function append_residues!(
    chain::Chain,
    dihedral_angles::AbstractVector{<:Real},
    bond_lengths::AbstractVector{<:Real} = IDEALIZED_BOND_LENGTHS,
    bond_angles::AbstractVector{<:Real} = IDEALIZED_BOND_ANGLES,
)
end=#