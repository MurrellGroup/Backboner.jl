export assign_secondary_structure!, assign_secondary_structure

import AssigningSecondaryStructure: assign_secondary_structure!, assign_secondary_structure

"""
    assign_secondary_structure!(protein)

Uses a simplified version of DSSP to fill the secondary structure vector of each chain with '-' (coil/loop), 'H' (helix), and 'E' (strand).
"""
function assign_secondary_structure!(protein::AbstractVector{ProteinChain})
    ss_vectors = assign_secondary_structure(NCaCO_coords.(protein))
    for (chain, ssvector) in zip(protein, ss_vectors)
        @assert length(chain.ssvector) == length(ssvector)
        chain.ssvector .= ssvector
    end
    return protein
end

"""
    assign_secondary_structure(protein)

Returns a new protein with secondary structure assigned.
"""
function assign_secondary_structure(protein::AbstractVector{ProteinChain})
    new_protein = deepcopy(protein)
    assign_secondary_structure!(new_protein)
    return new_protein
end

export has_assigned_ss

has_assigned_ss(ssvector::Vector{Char}) = all(!=(' '), ssvector)
has_assigned_ss(chain::ProteinChain) = has_assigned_ss(chain.ssvector)
has_assigned_ss(protein::AbstractVector{ProteinChain}) = all(has_assigned_ss, protein)