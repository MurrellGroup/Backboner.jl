const BACKBONE_BOND_LENGTHS = Float32[1.45775, 1.52307, 1.33208]
const BACKBONE_BOND_ANGLES = Float32[1.93731, 2.03926, 2.12710]

function Backboner.idealize(chain::Protein.Chain, args...; kwargs...)
    return Protein.Chain(
        chain.id,
        Backboner.idealize(chain.backbone, BACKBONE_BOND_LENGTHS, BACKBONE_BOND_ANGLES, args...; kwargs...),
        modelnum=chain.modelnum,
        resnums=chain.resnums,
        aavector=chain.aavector,
        ssvector=chain.ssvector,
    )
end