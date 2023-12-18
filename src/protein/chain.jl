export ProteinChain

"""
    ProteinChain <: AbstractVector{Residue}

A `ProteinChain` represents a chain of a protein, and is a vector of `Residue`s, which are instantiated from indexing the chain.

## Fields
- `id::AbstractString`: A string identifier (usually a single letter).
- `backbone::Backbone{3}`: An backbone with 3 atoms per residue (N, Ca, C), storing the coordinates of backbone atoms in a 3x3xL array.
- `aavector::Vector{Char}`: storing the amino acid sequence.
- `ssvector::Vector{Char}`: storing the secondary structure.
"""
struct ProteinChain <: AbstractVector{Residue}
    id::AbstractString
    backbone::Backbone{3}
    aavector::Vector{Char}
    ssvector::Vector{Char}

    function ProteinChain(
        id::AbstractString,
        backbone::Backbone{3};
        aavector::Vector{Char} = fill('G', size(backbone, 3)),
        ssvector::Union{Vector{Char}, Vector{<:Integer}} = fill(' ', size(backbone, 3)),
    )
        @assert size(backbone, 3) == length(aavector) == length(ssvector) "backbone, aavector, and ssvector must have the same length"
        ssvector isa Vector{<:Integer} && (ssvector = get.(('-', 'H', 'E'), ssvector, ' '))

        return new(id, backbone, aavector, ssvector)
    end

    ProteinChain(backbone::Backbone{3}; kwargs...) = ProteinChain("_", backbone; kwargs...) 
end

@inline Base.:(==)(chain1::ProteinChain, chain2::ProteinChain) = chain1.id == chain2.id && chain1.backbone == chain2.backbone && chain1.ssvector == chain2.ssvector
@inline Base.length(chain::ProteinChain) = size(chain.backbone, 3)
@inline Base.size(chain::ProteinChain) = Tuple(length(chain))
@inline Base.getindex(chain::ProteinChain, i::Integer) = Residue(i, chain.backbone, chain.aavector[i], chain.ssvector[i])

Base.summary(chain::ProteinChain) = "ProteinChain $(chain.id) with $(length(chain)) residue$(length(chain) == 1 ? "" : "s")"
Base.show(io::IO, chain::ProteinChain) = print(io, summary(chain))

@inline Base.getindex(protein::AbstractVector{ProteinChain}, i::AbstractString) = protein[findfirst(c -> c.id == i, protein)]
@inline Base.getindex(protein::AbstractVector{ProteinChain}, i::Symbol) = protein[String(i)]

export nitrogen_alphacarbon_distances, alphacarbon_carbonyl_distances, carbonyl_nitrogen_distances

nitrogen_alphacarbon_distances(backbone::Backbone{3}) = get_atom_distances(backbone, 1, 2, 0)
alphacarbon_carbonyl_distances(backbone::Backbone{3}) = get_atom_distances(backbone, 2, 3, 0)
carbonyl_nitrogen_distances(backbone::Backbone{3}) = get_atom_distances(backbone, 3, 1, 1)

"""
    nitrogen_alphacarbon_distances(chain::ProteinChain)

Calculate the distances between all pairs of contiguous nitrogen and alpha-carbon atoms in a chain.
Returns a vector of distances of length `length(chain)`.
"""
nitrogen_alphacarbon_distances(chain::ProteinChain) = nitrogen_alphacarbon_distances(chain.backbone)

"""
    alphacarbon_carbonyl_distances(chain::ProteinChain)

Calculate the distances between all pairs of contiguous alpha-carbon and carbonyl atoms in a chain.
Returns a vector of distances of length `length(chain)`.
"""
alphacarbon_carbonyl_distances(chain::ProteinChain) = alphacarbon_carbonyl_distances(chain.backbone)

"""
    carbonyl_nitrogen_distances(chain::ProteinChain)

Calculate the distances between all pairs of contiguous carbonyl and nitrogen atoms in a chain.
Returns a vector of distances of length `length(chain) - 1`.
"""
carbonyl_nitrogen_distances(chain::ProteinChain) = carbonyl_nitrogen_distances(chain.backbone)