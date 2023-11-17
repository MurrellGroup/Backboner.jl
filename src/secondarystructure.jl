export SecondaryStructure, Unassigned, Loop, Helix, Strand, has_missing_ss

@enum SecondaryStructure Unassigned Loop Helix Strand

has_missing_ss(ssv::AbstractVector{SecondaryStructure}) = Unassigned in ssv