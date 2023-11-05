export AbstractChain

"""
    AbstractChain{A, T} <: AbstractArray{T, 3}

An abstract type for a chain of residues.
All chain types have a `coords` field that is a 3xAxL matrix of coordinates.
"""
abstract type AbstractChain{A, T <: Real} <: AbstractArray{T, 3} end

@inline Base.size(chain::AbstractChain) = size(chain.coords)
@inline Base.length(chain::AbstractChain) = size(chain, 3)
@inline Base.getindex(chain::AbstractChain, i, j, k) = chain.coords[i,j,k]
@inline Base.getindex(chain::AbstractChain, i::Integer) = view(chain.coords, :, :, i)

Base.summary(chain::AbstractChain{A, T}) where {A, T} = "Chain $(chain.id) with $(length(chain)) residues"
Base.show(io::IO, chain::AbstractChain{A, T}) where {A, T} = print(io, summary(chain))

include("coordmatrix.jl")
include("chain.jl")
include("segment.jl")