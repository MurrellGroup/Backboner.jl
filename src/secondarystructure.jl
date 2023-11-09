export SecondaryStructure, MiSSing, Loop, Helix, Strand, has_missing_ss

@enum SecondaryStructure MiSSing Loop Helix Strand

has_missing_ss(ssv::AbstractVector{SecondaryStructure}) = MiSSing in ssv