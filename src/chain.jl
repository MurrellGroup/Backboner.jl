export Chain

"""
    Chain <: AbstractVector{Residue}

A chain has an identifier (usually a single letter) and holds the backbone atom coordinates, amino acid sequence, and secondary structures of a protein chain. 
"""
struct Chain <: AbstractVector{Residue}
    id::AbstractString
    backbone::Backbone{4}
    aavector::Vector{Char}
    ssvector::Vector{Char}

    function Chain(
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

    Chain(backbone::Backbone; kwargs...) = Chain("_", backbone; kwargs...) 
end

@inline Base.:(==)(chain1::Chain, chain2::Chain) = chain1.id == chain2.id && chain1.backbone == chain2.backbone && chain1.ssvector == chain2.ssvector
@inline Base.length(chain::Chain) = length(chain.backbone)
@inline Base.size(chain::Chain) = (length(chain),)
@inline Base.getindex(chain::Chain, i::Integer) = Residue(i, chain.backbone, chain.aavector[i], chain.ssvector[i])

Base.summary(chain::Chain) = "Chain $(chain.id) with $(length(chain)) residue$(length(chain) == 1 ? "" : "s")"
Base.show(io::IO, chain::Chain) = print(io, summary(chain))

has_assigned_ss(chain::Chain) = has_assigned_ss(chain.ssvector)
