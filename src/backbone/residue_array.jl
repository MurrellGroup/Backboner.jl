# (residues, atoms, 3)
struct ResidueArray{A, T <: Real} <: AbstractArray{T, 3}
    coords::AbstractArray{T, 3}

    function ResidueArray(coords::AbstractArray{T, 3}) where T <: Real
        A = size(coords, 2)
        @assert size(coords, 3) == 3 "coords must have 3 coordinates per atom"
        return new{A, T}(coords)
    end

    function ResidueArray{A}(coords::AbstractArray{T, 3}) where {A, T <: Real}
        @assert size(coords, 2) == A "coords must have $A columns"
        @assert size(coords, 3) == 3 "coords must have 3 coordinates per atom"
        return new{A, T}(coords)
    end
end

@inline Base.size(ra::ResidueArray) = size(ra.coords)
@inline Base.getindex(ra::ResidueArray, I::Vararg{Integer}) = ra.coords[I...]
@inline Base.getindex(ra::ResidueArray, r::UnitRange{<:Integer}) = ResidueArray(view(ra.coords, r, :, :))
@inline nresidues(ra::ResidueArray) = size(ra, 1)

@inline function _ith_column_coords(ra::ResidueArray{A}, i::Int) where A
    @assert i <= A "i must be less than or equal to the number of atoms per residue"
    L = nresidues(ra)
    return view(reshape(ra.coords, :, A), (i-1)*L+1:i*L, :)' # (3, L)
end