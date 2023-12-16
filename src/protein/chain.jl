export ProteinChain

"""
    ProteinChain <: AbstractVector{Residue}

A chain has an identifier (usually a single letter) and holds the backbone atom coordinates, amino acid sequence, and secondary structures of a protein chain. 
"""
struct ProteinChain <: AbstractVector{Residue}
    id::AbstractString
    backbone::Backbone{3}
    aavector::Vector{Char}
    ssvector::Vector{Char}

    function ProteinChain(
        id::AbstractString,
        backbone::Backbone{3};
        aavector::Vector{Char} = fill('G', length(backbone)),
        ssvector::Union{Vector{Char}, Vector{<:Integer}} = fill(' ', length(backbone)),
    )
        @assert length(backbone) == length(aavector) == length(ssvector) "backbone, aavector, and ssvector must have the same length"
        ssvector isa Vector{<:Integer} && (ssvector = get.(('-', 'H', 'E'), ssvector, ' '))

        return new(id, backbone, aavector, ssvector)
    end

    ProteinChain(backbone::Backbone{3}; kwargs...) = ProteinChain("_", backbone; kwargs...) 
end

@inline Base.:(==)(chain1::ProteinChain, chain2::ProteinChain) = chain1.id == chain2.id && chain1.backbone == chain2.backbone && chain1.ssvector == chain2.ssvector
@inline Base.length(chain::ProteinChain) = length(chain.backbone)
@inline Base.size(chain::ProteinChain) = (length(chain),)
@inline Base.getindex(chain::ProteinChain, i::Integer) = Residue(i, chain.backbone, chain.aavector[i], chain.ssvector[i])

Base.summary(chain::ProteinChain) = "ProteinChain $(chain.id) with $(length(chain)) residue$(length(chain) == 1 ? "" : "s")"
Base.show(io::IO, chain::ProteinChain) = print(io, summary(chain))

@inline Base.getindex(protein::AbstractVector{ProteinChain}, id::AbstractString) = protein[findfirst(c -> c.id == id, protein)]

"""
    nitrogen_alphacarbon_distances(backbone::Backbone)

Calculate the distances between all pairs of contiguous nitrogen and alpha-carbon atoms in a backbone.
Returns a vector of distances of length `length(backbone)`.
"""
nitrogen_alphacarbon_distances(backbone::Backbone{3}) = atom_distances(backbone, 1, 2)

"""
    alphacarbon_carbonyl_distances(backbone::Backbone)

Calculate the distances between all pairs of contiguous alpha-carbon and carbonyl atoms in a backbone.
Returns a vector of distances of length `length(backbone)`.
"""
alphacarbon_carbonyl_distances(backbone::Backbone{3}) = atom_distances(backbone, 2, 3)

"""
    carbonyl_nitrogen_distances(backbone::Backbone)

Calculate the distances between all pairs of contiguous carbonyl and nitrogen atoms in a backbone.
Returns a vector of distances of length `length(backbone) - 1`.
"""
carbonyl_nitrogen_distances(backbone::Backbone{3}) = atom_distances(backbone, 3, 1, 1)