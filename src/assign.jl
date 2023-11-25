export assign_secondary_structure!, assign_secondary_structure

import AssigningSecondaryStructure: assign_secondary_structure!, assign_secondary_structure

"""
    assign_secondary_structure!(protein)

Uses a simplified version of DSSP to fill the secondary structure vector of each chain with Loop, Helix, and Strand.
"""
function assign_secondary_structure!(protein::Protein)
    ss_num_vectors = assign_secondary_structure([chain.backbone.coords for chain in protein])
    for (chain, ss_num_vector) in zip(protein, ss_num_vectors)
        ssvec = SecondaryStructure.(ss_num_vector)
        @assert length(chain.ssvec) == length(ssvec)
        chain.ssvec .= ssvec
    end
    return protein
end

"""
    assign_secondary_structure(protein)

Returns a new protein with secondary structure assigned.
"""
function assign_secondary_structure(protein::Protein)
    new_protein = deepcopy(protein)
    assign_secondary_structure!(new_protein)
    return new_protein
end