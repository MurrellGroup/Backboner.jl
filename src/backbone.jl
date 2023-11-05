export Backbone, remove_column, has_missing_ss

"""
    Backbone{A, T} <: AbstractVector{Chain{A, T}}

A wrapper for a vector of chains.
Chains can be accessed by index or by ID.
"""
struct Backbone{A, T <: Real} <: AbstractVector{Chain{A, T}}
    chains::Vector{Chain{A, T}}
    id_dict::Dict{String, Chain{A, T}}

    function Backbone(chains::Vector{Chain{A, T}}) where {A, T}
        @assert length(unique([chain.id for chain in chains])) == length(chains)
        id_dict = Dict{String, Chain{A, T}}(chain.id => chain for chain in chains)
        return new{A, T}(chains, id_dict)
    end
end

@inline Base.size(bb::Backbone) = size(bb.chains)
@inline Base.length(bb::Backbone) = length(bb.chains)
@inline Base.getindex(bb::Backbone, i) = bb.chains[i]
@inline Base.getindex(bb::Backbone, id::AbstractString) = bb.id_dict[String(id)]

Base.summary(bb::Backbone{A, T}) where {A, T} = "Backbone{$A, $T} with $(length(bb)) chains"

has_missing_ss(backbone::Backbone) = any(has_missing_ss, backbone.chains)

function remove_column(backbone::Backbone{A, T}, i::Integer) where {A, T}
    new_chains = Vector{Chain{A-1, T}}()
    for chain in backbone
        push!(new_chains, remove_column(chain, i))
    end
    return Backbone(new_chains)
end