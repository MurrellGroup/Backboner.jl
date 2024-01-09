# module PDB

import BioStructures

const STANDARD_RESIDUE_ANGSTROM = [
    -1.066  -0.200   1.266;
     0.645  -0.527  -0.118;
     0.000   0.000   0.000;
] #  N       Ca      C

const AMINOACIDS = Set("ACDEFGHIKLMNPQRSTVWY")
const BACKBONE_ATOM_NAMES = ["N", "CA", "C"]

oneletter_resname(threeletter::AbstractString) = Char(get(BioStructures.threeletter_to_aa, threeletter, 'X'))

oneletter_resname(res::BioStructures.AbstractResidue) = oneletter_resname(BioStructures.resname(res))

backboneselector(at::BioStructures.AbstractAtom) = BioStructures.atomnameselector(at, BACKBONE_ATOM_NAMES)

proteindesignselector(res::BioStructures.AbstractResidue) =
    oneletter_resname(res) ∈ AMINOACIDS &&
    BioStructures.countatoms(res, backboneselector)==3 && 
    BioStructures.standardselector(res) &&
    !BioStructures.disorderselector(res)

atomcoords(res, ATOM_NAME)::Vector = BioStructures.collectatoms(res, at -> BioStructures.atomnameselector(at, [ATOM_NAME])) |> only |> BioStructures.coords

backbonecoords(res::BioStructures.AbstractResidue)::Matrix = stack(AT -> atomcoords(res, AT), BACKBONE_ATOM_NAMES)

getresidues(chain::BioStructures.Chain, selectors... = proteindesignselector) = BioStructures.collectresidues(chain, selectors...)

issmall(chain::BioStructures.Chain; threshold = 5) = length(getresidues(chain)) <= threshold

Backboner.Backbone(chain::BioStructures.Chain) = Backbone(convert(Matrix{Float32}, mapreduce(backbonecoords, hcat, getresidues(chain))))

function Backboner.Backbone(struc::BioStructures.ProteinStructure)
    backbones = Backbone{Float32}[]
    for model in struc, chain in model
        issmall(chain) || push!(backbones, Backbone(chain))
    end

    return Backbone(mapreduce(b -> b.coords, hcat, backbones, init=Matrix{Float32}(undef, 3, 0)))
end

Backboner.Frames(struc::BioStructures.ProteinStructure) = Frames(Backbone(struc), STANDARD_RESIDUE_ANGSTROM)

sequence(chain::BioStructures.Chain) = join(map(oneletter_resname, getresidues(chain)))

function sequence(struc::BioStructures.ProteinStructure)
    _sequence = String[]
    for model in struc, chain in model
        issmall(chain) || push!(_sequence, sequence(chain))
    end
    return join(_sequence)
end

modelspecific_chainid(chain::BioStructures.Chain) = findfirst(==(BioStructures.chainid(chain)), BioStructures.chainids(BioStructures.model(chain)))

modelnumber(chain::BioStructures.Chain) = BioStructures.modelnumber(chain)

function allchainids(struc::BioStructures.ProteinStructure)
    _allchainids = Tuple{Int, Int}[]
    for model in struc, chain in model
        issmall(chain) || push!(_allchainids, (modelspecific_chainid(chain), modelnumber(chain)))
    end
    return sort(_allchainids)
end

chainid(chain::BioStructures.Chain, all_chainids::Vector) = Int32(searchsortedfirst(all_chainids, (modelspecific_chainid(chain), modelnumber(chain))))

chainids(chain::BioStructures.Chain, all_chainids) = fill(chainid(chain, all_chainids), length(getresidues(chain)))

function chainids(struc::BioStructures.ProteinStructure)
    _chainids = Int32[]

    all_chainids = allchainids(struc)
    for model in struc, chain in model
        issmall(chain) || append!(_chainids, chainids(chain, all_chainids))
    end
    return _chainids
end

"""
SEQRES   1 A  548  MET ALA GLN LEU SER GLY GLN PRO VAL VAL ILE LEU PRO          
SEQRES   2 A  548  GLU GLY THR GLN ARG TYR VAL GLY ARG ASP ALA GLN ARG          
SEQRES   3 A  548  LEU ASN ILE LEU ALA ALA ARG ILE ILE ALA GLU THR VAL          
SEQRES   4 A  548  ARG THR THR LEU GLY PRO LYS GLY MET ASP LYS MET LEU          
SEQRES   5 A  548  VAL ASP SER LEU GLY ASP ILE VAL VAL THR ASN ASP CYS          
SEQRES   6 A  548  ALA THR ILE LEU ASP LYS ILE ASP LEU GLN HIS PRO ALA          
SEQRES   7 A  548  ALA LYS MET MET VAL GLU VAL ALA LYS THR GLN ASP LYS 
"""

function read_seqres_line(line::String)
    cols = split(line)
    _chainid = cols[3]
    residues = oneletter_resname.(cols[5:end]) # a gap will become an X

    return _chainid, residues
end

function read_chainwise_seqres(pdbfile::String)
    chainwise_seqres = Dict{String, Vector{Char}}()

    for line in eachline(pdbfile)
        if startswith(line, "SEQRES")
            _chainid, residues = read_seqres_line(line)
            append!(get!(chainwise_seqres, _chainid, Char[]), residues)
        end
    end

    return Dict{String, String}(map(k -> k => join(chainwise_seqres[k]), collect(keys(chainwise_seqres))))
end

function resnumbers(chain::BioStructures.Chain, chainwise_seqres)
    backbone = Backboner.Backbone(chain)
    resnumbers(sequence(chain), chainwise_seqres[BioStructures.chainid(chain)], backbone)
end

function resnumbers(struc::BioStructures.ProteinStructure, chainwise_seqres)
    _resnumbers = Int[]
    for model in struc, chain in model
        issmall(chain) || append!(_resnumbers, resnumbers(chain, chainwise_seqres))
    end
    return _resnumbers
end
# end