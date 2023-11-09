export Segment, extend_segment, segments

"""
    Segment{SS, A, T} <: AbstractChain{A, T}

A segment of a chain with the same secondary structure.
"""
struct Segment{SS, A, T} <: AbstractChain{A, T}
    chain::Chain{A, T}
    range::UnitRange{Int}
    coords::AbstractArray{T, 3}

    function Segment{SS}(chain::Chain{A, T}, r::UnitRange{Int}) where {SS, A, T}
        @assert SS == MiSSing || all(==(SS), view(chain.ssvector, r)) "All residues in the '$SS' segment must have the '$SS' secondary structure"
        coords = view(chain.coords, :, :, r)
        return new{SS, A, T}(chain, r, coords)
    end
end

@inline Base.getindex(segment::Segment{SS}, r::UnitRange{Int}) where SS = Segment{SS}(segment.chain, segment.range[r])

"""
    extend_segment(segment, range)

Returns a segment of the parent chain, extended to the given range.
If `segment` covers indices 3:4 of the parent chain, then `extend_segment(segment, 0:3)` will return a segment that covers indices 2:5, since 0 is one less than 1, and 3 is one more than 4.
`extend_segment(segment, 1:2)` would therefore return the same segment as `segment`.
This function is useful if one wishes to access the coordinates of the atoms of the parent chain that are adjacent to the segment.
!!! note
    The new segment will have missing secondary structure, since segments are defined by the secondary structure of the residues.
"""
@inline function extend_segment(segment::Segment{SS}, r::UnitRange{Int}) where SS
    offset = segment.range.start - 1
    parent_vec_range = r .+ offset
    checkbounds(segment.chain, parent_vec_range)
    return Segment{MiSSing}(segment.chain, parent_vec_range)
end

Base.summary(segment::Segment{SS}) where SS = "$SS Segment of Chain $(segment.chain.id) with $(length(segment)) residues"

"""
    segments(chain)

Returns an array of segments of a chain.
The segments are defined by the secondary structure of the residues.
A chain with missing secondary structure information will throw an error.
"""
function segments(chain::Chain)
    has_missing_ss(chain) && error("Chain $(chain.id) has missing secondary structure information")
    ssvector = chain.ssvector
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