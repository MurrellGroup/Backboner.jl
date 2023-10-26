struct Segment{C}
    range::UnitRange{Int}
end

function segmenter(ss::AbstractVector{ASS.SSClass})
    start_idx = 1
    end_idx = 1
    segments = Segment[]
    for (i, class) in enumerate(ss)
        if class != ss[start_idx]
            push!(segments, Segment{ss[start_idx]}(start_idx:end_idx))
            start_idx = i
        end
        end_idx = i
    end
    push!(segments, Segment{ss[start_idx]}(start_idx:end_idx))
    return segments
end

export segmenter