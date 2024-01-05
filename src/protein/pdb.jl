export readpdb

import BioStructures

const AMINOACIDS = Set("ACDEFGHIKLMNPQRSTVWY")
const BACKBONE_ATOM_NAMES = ["N", "CA", "C"]

oneletter_resname(res::BioStructures.Residue) = Char(get(BioStructures.threeletter_to_aa, BioStructures.resname(res), '-'))

backboneselector(at::BioStructures.AbstractAtom) = BioStructures.atomnameselector(at, BACKBONE_ATOM_NAMES)

proteindesignselector(res::BioStructures.Residue) =
    oneletter_resname(res) ∈ AMINOACIDS &&
    BioStructures.countatoms(res, backboneselector)==3 && 
    BioStructures.standardselector(res) &&
    !BioStructures.disorderselector(res)

atomcoords(res, ATOM_NAME)::Vector = BioStructures.collectatoms(res, at -> BioStructures.atomnameselector(at, [ATOM_NAME])) |> only |> BioStructures.coords

backbonecoords(res::BioStructures.Residue)::Matrix = stack(AT -> atomcoords(res, AT), BACKBONE_ATOM_NAMES)

getresidues(chain::BioStructures.Chain) = BioStructures.collectresidues(chain, proteindesignselector)

Base.isempty(chain::BioStructures.Chain) = isempty(getresidues(chain))

Backboner.Backbone(chain::BioStructures.Chain) = Backbone(mapreduce(backbonecoords, hcat, getresidues(chain)))

function Backboner.Backbone(struc::BioStructures.ProteinStructure)
    backbones = Backbone[]
    for model in struc, chain in model
        isempty(chain) || push!(backbones, Backbone(chain))
    end
    if isempty(backbones)
        return Backbone{Float64}(undef, 0)
    else
        return Backbone(mapreduce(b -> b.coords, hcat, backbones))
    end
end

sequence(chain::BioStructures.Chain) = join(map(oneletter_resname, getresidues(chain)))

function sequence(struc::BioStructures.ProteinStructure)
    _sequence = String[]
    for model in struc, chain in model
        isempty(chain) || push!(_sequence, sequence(chain))
    end
    return join(_sequence)
end

function modelspecific_chainid(chain::BioStructures.Chain)
    str = BioStructures.chainid(chain)
    if length(str) <= 4
        throw(ArgumentError("ChainID is too long for structure $(BioStructures.structurename(chain)), chain $(str)"))
    end
    int_bits = join(map(bitstring, collect(codeunits(str))))
    return parse(Int32, int_bits)
end

modelnumber(chain::BioStructures.Chain) = Int32(BioStructures.modelnumber(chain))

chainid(chain::BioStructures.Chain) = parse(Int64, join([
        modelnumber(chain)::Int32 |> bitstring,
        modelspecific_chainid(chain)::Int32 |> bitstring
    ]))

chainids(chain::BioStructures.Chain) = fill(chainid(chain), length(getresidues(chain)))

function chainids(struc::BioStructures.ProteinStructure)
    _chainids = Vector{Int64}[]
    for model in struc, chain in model
        isempty(chain) || push!(_chainids, chainids(chain))
    end
    return vcat(_chainids)
end