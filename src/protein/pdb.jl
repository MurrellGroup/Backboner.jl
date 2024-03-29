export readpdb, writepdb

import BioStructures

const AMINOACIDS = Set("ACDEFGHIKLMNPQRSTVWY")
const BACKBONE_ATOM_NAMES = ["N", "CA", "C"]

oneletter_resname(threeletter::AbstractString) = Char(get(BioStructures.threeletter_to_aa, threeletter, 'X'))
oneletter_resname(res::BioStructures.AbstractResidue) = oneletter_resname(BioStructures.resname(res))

backboneselector(at::BioStructures.AbstractAtom) = BioStructures.atomnameselector(at, BACKBONE_ATOM_NAMES)
backboneselector(res::BioStructures.AbstractResidue) =
    oneletter_resname(res) in AMINOACIDS &&
    BioStructures.countatoms(res, backboneselector) == 3 &&
    BioStructures.standardselector(res) &&
    !BioStructures.disorderselector(res)

getresidues(chain::BioStructures.Chain, res_selector) = BioStructures.collectresidues(chain, res_selector)

atomcoords(res, atom_name)::Vector = BioStructures.collectatoms(res, atom -> BioStructures.atomnameselector(atom, [atom_name])) |> only |> BioStructures.coords

backbonecoords(res::BioStructures.AbstractResidue)::Matrix = stack(AT -> atomcoords(res, AT), BACKBONE_ATOM_NAMES)

protein_backbone(residues::Vector{<:BioStructures.AbstractResidue}) = Backboner.Backbone(convert(Matrix{Float32}, mapreduce(backbonecoords, hcat, residues; init=Matrix{Float32}(undef, 3, 0))))
protein_backbone(chain::BioStructures.Chain, res_selector) = Backboner.Backbone(getresidues(chain, res_selector))

aminoacid_sequence(residues::Vector{<:BioStructures.AbstractResidue}) = map(oneletter_resname, residues)
aminoacid_sequence(chain::BioStructures.Chain, res_selector) = aminoacid_sequence(getresidues(chain, res_selector))

function Protein.Chain(residues::Vector{<:BioStructures.AbstractResidue})
    id = only(unique(BioStructures.chainid.(residues)))
    backbone = protein_backbone(residues)
    aavector = aminoacid_sequence(residues)
    resnums = BioStructures.resnumber.(residues)
    modelnum = only(unique(BioStructures.modelnumber.(residues)))
    return Protein.Chain(id, backbone; resnums, aavector, modelnum)
end

function Protein.Chain(chain::BioStructures.Chain; res_selector=backboneselector)
    residues = getresidues(chain, res_selector)
    isempty(residues) && return Protein.Chain(
        BioStructures.chainid(chain), Backbone{Float32}(undef, 0); resnums=Int[], aavector=Char[], modelnum=BioStructures.modelnumber(chain))
    return Protein.Chain(residues)
end

function chains(struc::BioStructures.ProteinStructure; res_selector=backboneselector)
    chains = Protein.Chain[]
    for model in struc, chain in model
        isempty(chain) || push!(chains, Protein.Chain(chain, res_selector=res_selector))
    end
    return chains
end

"""
    readpdb(pdbfile::String)

Loads a protein (represented as a `Vector{Protein.Chain}`) from a PDB file.
Assumes that each residue starts with three atoms: N, CA, C.
"""
function readpdb(pdbfile::String)
    struc = read(pdbfile, BioStructures.PDB)
    return chains(struc)
end


import PDBTools

"""
    writepdb(protein::Vector{Protein.Chain}, filename)

Write a protein (represented as a `Vector{Protein.Chain}`s) to a PDB file.
"""
function writepdb(protein::Vector{Chain}, filename, header=:auto, footer=:auto)
    atoms = PDBTools.Atom[]
    index = 0
    residue_index = 0
    for chain in protein
        coords = ncaco_coords(chain)
        for (residue, res_coords) in zip(chain, eachslice(coords, dims=3))
            resname = get(threeletter_aa_names, residue.aa, "XXX")
            residue_index += 1
            for (name, atomcoords) in zip(["N", "CA", "C", "O"], eachcol(res_coords))
                index += 1
                atom = PDBTools.Atom(
                    index = index,
                    name = name,
                    resname = resname,
                    chain = chain.id,
                    resnum = residue.num,
                    residue = residue_index,
                    x = atomcoords[1],
                    y = atomcoords[2],
                    z = atomcoords[3],
                )
                push!(atoms, atom)
            end
        end
    end
    PDBTools.writePDB(atoms, filename, header=header, footer=footer)
end
