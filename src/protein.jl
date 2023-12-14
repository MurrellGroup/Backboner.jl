@inline Base.getindex(protein::AbstractVector{ProteinChain}, id::AbstractString) = protein[findfirst(c -> c.id == id, protein)]

has_assigned_ss(protein::AbstractVector{ProteinChain}) = all(has_assigned_ss, protein)