export AbstractChain, atom_coord_matrix, unroll_atoms

abstract type AbstractChain{A, T <: Real} <: AbstractArray{T, 3} end

@inline Base.size(chain::AbstractChain) = size(chain.coords)
@inline Base.length(chain::AbstractChain) = size(chain, 3)
@inline Base.getindex(chain::AbstractChain, i, j, k) = chain.coords[i,j,k]
@inline Base.getindex(chain::AbstractChain, i::Integer) = view(chain.coords, :, :, i)

Base.summary(chain::AbstractChain{A, T}) where {A, T} = "Chain $(chain.id) with $(length(chain)) residues"
Base.show(io::IO, chain::AbstractChain{A, T}) where {A, T} = print(io, summary(chain))

"""
    atom_coord_matrix(chain, i)

Returns the coordinates of specific columns of atoms in a chain.
"""
function atom_coord_matrix(chain::AbstractChain, i)
    return view(chain.coords, :, i, :)
end

"""
    unroll_atoms(chain, i)

Returns the coordinates of specific columns of atoms in a chain,
but unrolled into a 3xN matrix where N is the number of residues times
the number of columns selected (atoms selected per residue).
"""
function unroll_atoms(chain::AbstractChain, i=:)
    return reshape(atom_coord_matrix(chain, i), 3, :)
end

include("chain.jl")
include("segment.jl")