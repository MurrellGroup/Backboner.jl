#=--- Frame idealization ---=#

const STANDARD_RESIDUE_ANGSTROM = [
    -1.066  -0.200   1.266;
     0.645  -0.527  -0.118;
     0.000   0.000   0.000;
] #  N       Ca      C

#=--- Bond idealization ---=#

const BACKBONE_BOND_LENGTHS = Float32[1.45775, 1.52307, 1.33208]
const BACKBONE_BOND_ANGLES = Float32[1.93731, 2.03926, 2.12710]

function Backboner.idealize(chain::Protein.Chain, ideal_lengths=BACKBONE_BOND_LENGTHS, ideal_angles=BACKBONE_BOND_ANGLES; kwargs...)
    return Protein.Chain(
        chain.id,
        Backboner.idealize(chain.backbone, ideal_lengths, ideal_angles; kwargs...),
        modelnum=chain.modelnum,
        resnums=chain.resnums,
        aavector=chain.aavector,
        ssvector=chain.ssvector,
    )
end