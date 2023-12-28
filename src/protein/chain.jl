export Chain

"""
    Chain <: AbstractVector{Residue}

A `Chain` represents a chain of a protein, and is a vector of `Residue`s, which are instantiated from indexing the chain.

## Fields
- `id::AbstractString`: A string identifier (usually a single letter).
- `backbone::Backbone`: A backbone with a length divisible by 3, to ensure 3 atoms per residue (N, Ca, C).
- `aavector::Vector{Char}`: storing the amino acid sequence.
- `ssvector::Vector{Char}`: storing the secondary structure.
"""
struct Chain <: AbstractVector{Residue}
    id::AbstractString
    backbone::Backbone
    aavector::Vector{Char}
    ssvector::Vector{Char}

    function Chain(
        id::AbstractString,
        backbone::Backbone;
        aavector::Vector{Char} = fill('G', length(backbone) ÷ 3),
        ssvector::Union{Vector{Char}, Vector{<:Integer}} = fill(' ', length(backbone) ÷ 3),
    )
        @assert length(backbone) % 3 == 0
        @assert length(backbone) ÷ 3 == length(aavector) == length(ssvector) "backbone, aavector, and ssvector must have the same length"
        ssvector isa Vector{<:Integer} && (ssvector = get.(('-', 'H', 'E'), ssvector, ' '))

        return new(id, backbone, aavector, ssvector)
    end

    Chain(backbone::Backbone; kwargs...) = Chain("_", backbone; kwargs...) 
end

@inline Base.:(==)(chain1::Chain, chain2::Chain) = chain1.id == chain2.id && chain1.backbone == chain2.backbone && chain1.ssvector == chain2.ssvector
@inline Base.length(chain::Chain) = length(chain.backbone) ÷ 3
@inline Base.size(chain::Chain) = Tuple(length(chain))
@inline Base.getindex(chain::Chain, i::Integer) = Residue(i, chain.backbone, chain.aavector[i], chain.ssvector[i])

Base.summary(chain::Chain) = "Chain $(chain.id) with $(length(chain)) residue$(length(chain) == 1 ? "" : "s")"
Base.show(io::IO, chain::Chain) = print(io, summary(chain))

@inline Base.getindex(protein::AbstractVector{Chain}, id::AbstractString) = protein[findfirst(chain -> chain.id == id, protein)]
@inline Base.getindex(protein::AbstractVector{Chain}, id::Symbol) = protein[String(id)]

export nitrogens, alphacarbons, carbonyls

nitrogens(backbone::Backbone) = backbone[1:3:end]
alphacarbons(backbone::Backbone) = backbone[2:3:end]
carbonyls(backbone::Backbone) = backbone[3:3:end]

nitrogens(chain::Chain) = nitrogens(chain.backbone)
alphacarbons(chain::Chain) = alphacarbons(chain.backbone)
carbonyls(chain::Chain) = carbonyls(chain.backbone)
# oxygen_coords function in src/protein/oxygen.jl

export nitrogen_alphacarbon_distances, alphacarbon_carbonyl_distances, carbonyl_nitrogen_distances

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

phi_angles(bonds::ChainedBonds) = @view bonds.dihedrals[3:3:end]
psi_angles(bonds::ChainedBonds) = @view bonds.dihedrals[1:3:end]
omega_angles(bonds::ChainedBonds) = @view bonds.dihedrals[2:3:end]

# TODO: functions for getting beta-carbon and hydrogen atom positions? is that even possible without knowing AAs and molecular dynamics?

# FIXME: better idealized bond lengths and angles
# currently taken from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2810841/
#=const IDEALIZED_BOND_LENGTHS = [1.45, 1.52, 1.33]
const IDEALIZED_BOND_ANGLES = [121.7, 111.0, 117.2]

# TODO: append_residues! function. also get the Residue type sorted out. user shouldn't need to create it manually.
function append_residues!(
    chain::Chain,
    dihedral_angles::AbstractVector{<:Real},
    bond_lengths::AbstractVector{<:Real} = IDEALIZED_BOND_LENGTHS,
    bond_angles::AbstractVector{<:Real} = IDEALIZED_BOND_ANGLES,
)
end=#