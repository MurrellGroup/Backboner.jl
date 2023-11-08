export assign_secondary_structure!

import AssigningSecondaryStructure as ASS

"""
    assign_secondary_structure!(backbone)
"""
function assign_secondary_structure!(backbone::Backbone{4})
    ssvectors = ASS.assign_secondary_structure([chain.coords for chain in backbone])
    for (chain, ssvector) in zip(backbone, ssvectors)
        @assert length(chain.ssvector) == length(ssvector)
        chain.ssvector .= SecondaryStructure.(ssvector)
    end
end

function assign_secondary_structure!(backbone3::Backbone{3})
    backbone4 = backbone_with_oxygen(backbone3)
    assign_secondary_structure!(backbone4)
    for (chain3, chain4) in zip(backbone3, backbone4)
        chain3.ssvector .= chain4.ssvector
    end
end