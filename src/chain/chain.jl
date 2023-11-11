export Chain

# could be a wrapper for Vector{Residue} but that would overcomplicate things
"""
    Chain{T}

A chain has an identifier (usually a single letter) and holds the backbone atom coordinates and secondary structure of a protein chain. 
"""
struct Chain{T}
    id::AbstractString
    backbone::Backbone{4,T}
    ssvector::Vector{SecondaryStructure}

    function Chain(id::AbstractString, backbone::Backbone{4,T}, ssvector::Vector{SecondaryStructure}) where T
        @assert length(backbone) == length(ssvector) "backbone and ssvector must have the same length"
        return new{T}(id, backbone, ssvector)
    end

    function Chain(id::AbstractString, backbone::Backbone{4,T}) where T
        return Chain(id, backbone, fill(MiSSing, length(backbone)))
    end

    function Chain(id::AbstractString, backbone::Backbone{3})
        return Chain(id, add_oxygen_slice(backbone))
    end

    Chain(backbone::Backbone) = Chain("", backbone) 
end

@inline Base.:(==)(chain1::Chain, chain2::Chain) = chain1.id == chain2.id && chain1.backbone == chain2.backbone && (all(==(MiSSing), chain1.ssvector) || all(==(MiSSing), chain2.ssvector) || chain1.ssvector == chain2.ssvector)
@inline Base.length(chain::Chain) = length(chain.backbone)
@inline Base.size(chain::Chain) = (length(chain),)

Base.summary(chain::Chain) = "Chain $(chain.id) with $(length(chain)) residue$(length(chain) == 1 ? "" : "s")"
Base.show(io::IO, chain::Chain) = print(io, summary(chain))

has_missing_ss(chain::Chain) = has_missing_ss(chain.ssvector)

include("segment.jl")