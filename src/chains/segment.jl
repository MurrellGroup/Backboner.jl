export Segment, segments

"""
    Segment{SS, A, T} <: AbstractChain{A, T}

A segment of a chain with the same secondary structure.
"""
struct Segment{SS, A, T} <: AbstractChain{A, T}
    chain::Chain{A, T}
    range::UnitRange{Int}
    coords::AbstractArray{T, 3}

    function Segment{SS}(chain::Chain{A, T}, r::UnitRange{Int}) where {SS, A, T}
        @assert all(==(SS), view(chain.ssvector, r)) "All residues in the segment must have the same secondary structure"
        coords = view(chain.coords, :, :, r)
        return new{SS, A, T}(chain, r, coords)
    end
end

@inline Base.getindex(segment::Segment{SS}, r::UnitRange{Int}) where SS = Segment{SS}(segment.chain, segment.range[r])

Base.summary(segment::Segment{SS}) where SS = "$SS Segment of Chain $(segment.chain.id) with $(length(segment)) residues"

"""
    segments(chain)

Returns an array of segments of a chain.
The segments are defined by the secondary structure of the residues.
A chain with missing secondary structure information will throw an error.
"""
function segments(chain::Chain)
    ssvector = chain.ssvector
    any(==(MiSSing), ssvector) && error("Chain $(chain.id) has missing secondary structure information")
    start_idx = 1
    end_idx = 1
    segments = Segment[]
    for (i, ss) in enumerate(ssvector)
        if ss != ssvector[start_idx]
            push!(segments, Segment{ssvector[start_idx]}(chain, start_idx:end_idx))
            start_idx = i
        end
        end_idx = i
    end
    push!(segments, Segment{ssvector[start_idx]}(chain, start_idx:end_idx))
    return segments
end