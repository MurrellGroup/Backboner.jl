export assign_secondary_structure!

import AssigningSecondaryStructure as ASS

"""
    assign_secondary_structure!(protein)

Uses a simplified version of DSSP to fill the secondary structure vector of each chain with Loop, Helix, and Strand.
"""
function assign_secondary_structure!(protein::Protein)
    ss_num_vectors = ASS.assign_secondary_structure([chain.backbone.coords for chain in protein])
    for (chain, ss_num_vector) in zip(protein, ss_num_vectors)
        ssvec = SecondaryStructure.(ss_num_vector)
        @assert length(chain.ssvec) == length(ssvec)
        chain.ssvec .= ssvec
    end
end