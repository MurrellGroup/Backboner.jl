export Protein

"""
    Protein{T} <: AbstractVector{Chain{T}}

A wrapper for a vector of chains.
Chains can be accessed by index or by ID.
"""
struct Protein{T} <: AbstractVector{Chain{T}}
    chains::Vector{Chain{T}}
    id_dict::Dict{AbstractString, Chain{T}}

    function Protein(chains::Vector{Chain{T}}) where T
        @assert length(unique([chain.id for chain in chains])) == length(chains)
        id_dict = Dict{AbstractString, Chain{T}}(chain.id => chain for chain in chains)
        return new{T}(chains, id_dict)
    end
end

@inline Base.:(==)(protein1::Protein, protein2::Protein) = protein1.chains == protein2.chains
@inline Base.size(protein::Protein) = size(protein.chains)
@inline Base.length(protein::Protein) = length(protein.chains)
@inline Base.getindex(protein::Protein, i) = protein.chains[i]
@inline Base.getindex(protein::Protein, id::AbstractString) = protein.id_dict[String(id)]

Base.summary(protein::Protein{T}) where T = "Protein{$T} with $(length(protein)) chain$(length(protein) == 1 ? "" : "s")"

has_missing_ss(protein::Protein) = any(has_missing_ss, protein.chains)