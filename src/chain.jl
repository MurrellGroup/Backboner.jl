export ProteinChain

"""
    ProteinChain <: AbstractVector{Residue}

A chain has an identifier (usually a single letter) and holds the backbone atom coordinates, amino acid sequence, and secondary structures of a protein chain. 
"""
struct ProteinChain <: AbstractVector{Residue}
    id::AbstractString
    backbone::Backbone{4}
    aavector::Vector{Char}
    ssvector::Vector{Char}

    function ProteinChain(
        id::AbstractString,
        backbone::Backbone{N};
        aavector::Vector{Char} = fill('G', length(backbone)),
        ssvector::Union{Vector{Char}, Vector{<:Integer}} = fill(' ', length(backbone)),
    ) where N
        @assert N == 3 || N == 4 "backbone must have 3 or 4 atoms per residue"
        N == 3 && (backbone = add_oxygens(backbone))

        @assert length(backbone) == length(aavector) == length(ssvector) "backbone, aavector, and ssvector must have the same length"
        ssvector isa Vector{<:Integer} && (ssvector = get.("-HE", ssvector, ' '))

        return new(id, backbone, aavector, ssvector)
    end

    ProteinChain(backbone::Backbone; kwargs...) = ProteinChain("_", backbone; kwargs...) 
end

@inline Base.:(==)(chain1::ProteinChain, chain2::ProteinChain) = chain1.id == chain2.id && chain1.backbone == chain2.backbone && chain1.ssvector == chain2.ssvector
@inline Base.length(chain::ProteinChain) = length(chain.backbone)
@inline Base.size(chain::ProteinChain) = (length(chain),)
@inline Base.getindex(chain::ProteinChain, i::Integer) = Residue(i, chain.backbone, chain.aavector[i], chain.ssvector[i])

Base.summary(chain::ProteinChain) = "ProteinChain $(chain.id) with $(length(chain)) residue$(length(chain) == 1 ? "" : "s")"
Base.show(io::IO, chain::ProteinChain) = print(io, summary(chain))

has_assigned_ss(chain::ProteinChain) = has_assigned_ss(chain.ssvector)
