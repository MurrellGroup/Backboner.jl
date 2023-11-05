export Chain

"""
    Chain{A, T} <: AbstractChain{A, T}

A wrapper for a 3xAxL matrix of coordinates of atoms in a backbone chain.
"""
struct Chain{A, T} <: AbstractChain{A, T}
    id::AbstractString
    coords::AbstractArray{T, 3}
    ssvector::Vector{SecondaryStructure}

    function Chain(id::AbstractString, coords::AbstractArray{T, 3}, ss::Union{Nothing, Vector{SecondaryStructure}} = nothing) where T
        @assert size(coords, 1) == 3 "coords must have 3 coordinates per atom"
        A = size(coords, 2)
        ssvector = isnothing(ss) ? fill(MiSSing, size(coords, 3)) : ss
        return new{A, T}(id, coords, ssvector)
    end
end

@inline Base.getindex(chain::Chain, r::UnitRange{Int}) = Chain("$(chain.id)[$r]", view(chain.coords, :, :, r))

has_missing_ss(chain::Chain) = any(==(MiSSing), chain.ssvector)

function remove_column(chain::Chain{A, T}, i::Integer) where {A, T}
    @assert i <= A
    return Chain(chain.id, view(chain.coords, :, [1:i-1; i+1:A], :))
end