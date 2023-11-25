export Chain

# could be a wrapper for Vector{Residue} but that would overcomplicate things
"""
    Chain <: AbstractVector{Residue}

A chain has an identifier (usually a single letter) and holds the backbone atom coordinates, amino acid sequence, and secondary structures of a protein chain. 
"""
struct Chain <: AbstractVector{Residue}
    id::AbstractString
    backbone::Backbone{4}
    aaseq::Vector{Char}
    ssvec::Vector{SecondaryStructure}

    function Chain(
        id::AbstractString,
        backbone::Backbone{N};
        aaseq::Vector{Char} = fill('G', length(backbone)),
        ssvec::Union{Vector{SecondaryStructure}, Vector{<:Integer}} = fill(Unassigned, length(backbone)),
    ) where N
        @assert N == 3 || N == 4 "backbone must have 3 or 4 atoms per residue"
        N == 3 && (backbone = add_oxygens(backbone))

        @assert length(backbone) == length(aaseq) == length(ssvec) "backbone, aaseq, and ssvec must have the same length"
        ssvec isa Vector{<:Integer} && (ssvec = SecondaryStructure.(ssvec))

        return new(id, backbone, aaseq, ssvec)
    end

    Chain(backbone::Backbone; kwargs...) = Chain("_", backbone; kwargs...) 
end

@inline Base.:(==)(chain1::Chain, chain2::Chain) = chain1.id == chain2.id && chain1.backbone == chain2.backbone && chain1.ssvec == chain2.ssvec
@inline Base.length(chain::Chain) = length(chain.backbone)
@inline Base.size(chain::Chain) = (length(chain),)
@inline Base.getindex(chain::Chain, i::Integer) = Residue(i, chain.backbone, chain.aaseq[i], chain.ssvec[i])

Base.summary(chain::Chain) = "Chain $(chain.id) with $(length(chain)) residue$(length(chain) == 1 ? "" : "s")"
Base.show(io::IO, chain::Chain) = print(io, summary(chain))

has_missing_ss(chain::Chain) = has_missing_ss(chain.ssvec)
