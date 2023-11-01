export assign_secondary_structure!

import AssigningSecondaryStructure as ASS

function assign_secondary_structure!(backbone::Backbone{4})
    ssvectors = ASS.dssp([chain.coords for chain in backbone])
    for (chain, ssvector) in zip(backbone, ssvectors)
        @assert length(chain.ssvector) == length(ssvector)
        chain.ssvector .= SecondaryStructure.(ssvector)
    end
end