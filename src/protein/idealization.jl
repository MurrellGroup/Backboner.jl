### Frame idealization

const STANDARD_TRIANGLE_ANGSTROM = [
    -1.066  -0.200   1.266;
     0.645  -0.527  -0.118;
     0.000   0.000   0.000;
] #  N       Ca      C

# alias for backwards compatibility
const STANDARD_RESIDUE_ANGSTROM = STANDARD_TRIANGLE_ANGSTROM

### Bond idealization

const BACKBONE_BOND_LENGTHS = Float64[1.45775, 1.52307, 1.33208]
const BACKBONE_BOND_ANGLES = Float64[1.93731, 2.03926, 2.12710]

function Backboner.idealize(chain::Protein.Chain, ideal_lengths=BACKBONE_BOND_LENGTHS, ideal_angles=BACKBONE_BOND_ANGLES; kwargs...)
    return Protein.Chain(
        chain.id,
        Backboner.idealize(chain.backbone, ideal_lengths, ideal_angles; kwargs...),
        modelnum=chain.modelnum,
        resnums=chain.resnums,
        ins_codes=chain.ins_codes,
        aavector=chain.aavector,
        ssvector=chain.ssvector,
    )
end
