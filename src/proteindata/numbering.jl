#module Numbering

#using Backboner

"""
seqref = "ABCDEF"
subseq = "BCE"
resnumbers = [2, 3, 5]

seqref[resnumbers] == subseq || throw(ArgumentError("resnumbers do not match subseq - alignment is sketchy"))
"""

function resnumbers(subseq::String, refseq::String, backbone::Backbone;
    match = 5,
    mismatch = -5,
    gap_open = -10,
    gap_extend = -1,
    ideal_bond_length = 1.32,
    tol = 0.2
)
    contiguous = map(dist -> dist <= ideal_bond_length + tol, carbonyl_nitrogen_distances(backbone))

    any(contiguous) || throw(ArgumentError("No contiguous residues found"))

    MATCH = 1
    GAP = 2

    m = length(subseq)
    n = length(refseq)
    # Init dp
    # dp[i, j, MATCH] = highest score for matching subseq[1:i-1] with refseq[1:j-1], have not opened gap
    # dp[i, j, GAP] = highest score for having matched subseq[1:i-1] with refseq[1:j-1], have opened gap
    # We don't use half the dp array
    dp = fill(-Inf, m + 1, n + 1, 2)
    trace = fill(-1, m + 1, n + 1, 2)

    # If you match nothing with nothing, the score should be 0
    dp[1, :, MATCH] .= 0
    
    # Run dp
    for i in axes(dp, 1), j in axes(dp, 2)
        (i > 1 && j > 1) || continue
        match_cost = subseq[i - 1] == refseq[j - 1] ? match : mismatch

        if i == 2 || contiguous[i - 2]
            dp[i, j, MATCH] = match_cost + dp[i-1, j-1, MATCH]
            trace[i, j, MATCH] = MATCH
        else
            dp[i, j, MATCH] = match_cost + dp[i-1, j-1, GAP]
            trace[i, j, MATCH] = GAP
        end
        dp[i, j, GAP], trace[i, j, GAP] = findmax((dp[i, j-1, MATCH] + gap_open, dp[i, j-1, GAP] + gap_extend))
    end

    # Backtrack
    L = argmax(dp[end, 2:end, MATCH])
    
    mask = zeros(Bool, L)

    i = lastindex(dp, 1)
    j = 1 + L
    k = MATCH

    while i > 1
        prev_k = trace[i, j, k]
        j -= 1
        if k == MATCH
            i -= 1
            mask[j] = true
        end
        
        k = prev_k
    end

    resnumbers = findall(mask)
    
    return resnumbers
end

#=
get_contiguity(_carbonyl_nitrogen_distances; lo = 0.122f0, hi = 0.142f0) = map(x -> lo <= x <= hi, _carbonyl_nitrogen_distances)

function resnumbers(seq::String, refseq::String, backbone::AtomStack;
    match = 5,
    mismatch = -5,
    gap_open = -10,
    gap_extend = -1
)
    contiguous = backbone |> carbonyl_nitrogen_distances |> get_contiguity

    MATCH = 1
    GAP = 2

    m = length(seq)
    n = length(refseq)
    # Init dp
    # dp[i, j, MATCH] = highest score for matching seq[1:i-1] with refseq[1:j-1], have not opened gap
    # dp[i, j, GAP] = highest score for having matched seq[1:i-1] with refseq[1:j-1], have opened gap
    # We don't use half the dp array
    dp = fill(-Inf, m + 1, n + 1, 2)
    trace = fill(-1, m + 1, n + 1, 2)

    # If you match nothing with nothing, the score should be 0
    dp[1, :, MATCH] .= 0
    
    # Run dp
    for i in axes(dp, 1), j in axes(dp, 2)
        (i > 1 && j > 1) || continue
        match_cost = seq[i - 1] == refseq[j - 1] ? match : mismatch

        if i == 2 || contiguous[i - 2]
            dp[i, j, MATCH] = match_cost + dp[i-1, j-1, MATCH]
            trace[i, j, MATCH] = MATCH
        else
            dp[i, j, MATCH] = match_cost + dp[i-1, j-1, GAP]
            trace[i, j, MATCH] = GAP
        end
        dp[i, j, GAP], trace[i, j, GAP] = findmax((dp[i, j-1, MATCH] + gap_open, dp[i, j-1, GAP] + gap_extend))
    end

    # Backtrack
    L = argmax(dp[end, 2:end, MATCH])
    
    mask = zeros(Bool, L)

    i = lastindex(dp, 1)
    j = 1 + L
    k = MATCH

    while i > 1
        prev_k = trace[i, j, k]
        j -= 1
        if k == MATCH
            i -= 1
            mask[j] = true
        end
        
        k = prev_k
    end

    resnumbers = findall(mask)

    resnumbers_from_one = resnumbers .- minimum(resnumbers) .+ 1

    return resnumbers_from_one
end
=#

#end