export readcif, readpdb, writepdb

using BioStructures: PDBFormat, PDBXMLFormat, MMCIFFormat, MMTFFormat
const ProteinFileFormat = Union{PDBFormat, PDBXMLFormat, MMCIFFormat, MMTFFormat}
export PDBFormat, PDBXMLFormat, MMCIFFormat, MMTFFormat

# TODO: Vector{Chain} to MolecularStructure conversion to allow for writing to other formats than PDB
# temporary solution could be to write pdb, read+delete pdb, write e.g. cif

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

atomcoords(res, atom_name) = BioStructures.collectatoms(res, atom -> BioStructures.atomnameselector(atom, [atom_name])) |> only |> BioStructures.coords

backbonecoords(res::BioStructures.AbstractResidue)::Matrix = stack(AT -> atomcoords(res, AT), BACKBONE_ATOM_NAMES)

protein_backbone(residues::Vector{<:BioStructures.AbstractResidue}) = Backboner.Backbone(convert(Matrix{Float32}, mapreduce(backbonecoords, hcat, residues; init=Matrix{Float32}(undef, 3, 0))))
protein_backbone(chain::BioStructures.Chain, res_selector) = Backboner.Backbone(getresidues(chain, res_selector))

aminoacid_sequence(residues::Vector{<:BioStructures.AbstractResidue}) = map(oneletter_resname, residues)
aminoacid_sequence(chain::BioStructures.Chain, res_selector) = aminoacid_sequence(getresidues(chain, res_selector))

function get_residue_atoms(residues::Vector{<:BioStructures.AbstractResidue})
    residue_atoms = Vector{Atom}[]
    for res in residues
        push!(residue_atoms, Atom.(BioStructures.collectatoms(res, a -> BioStructures.standardselector(a) && !BioStructures.disorderselector(a))))
    end
    return residue_atoms
end

function Protein.Chain(residues::Vector{<:BioStructures.AbstractResidue})
    id = only(unique(BioStructures.chainid.(residues)))
    backbone = protein_backbone(residues)
    aavector = aminoacid_sequence(residues)
    resnums = BioStructures.resnumber.(residues)
    ins_codes = BioStructures.inscode.(residues)
    modelnum = only(unique(BioStructures.modelnumber.(residues)))
    residue_atoms = get_residue_atoms(residues)
    return Protein.Chain(backbone; id, resnums, ins_codes, aavector, modelnum, residue_atoms)
end

function Protein.Chain(chain::BioStructures.Chain; res_selector=backboneselector)
    residues = getresidues(chain, res_selector)
    isempty(residues) && return Protein.Chain(
        BioStructures.chainid(chain), Backbone(Matrix{Float32}(undef, 3, 0)); modelnum=BioStructures.modelnumber(chain))
    return Protein.Chain(residues)
end

function collectchains(struc::BioStructures.MolecularStructure; res_selector=backboneselector)
    chains = Protein.Chain[]
    for model in struc, chain in model
        isempty(chain) || push!(chains, Protein.Chain(chain, res_selector=res_selector))
    end
    return chains
end

"""
    readchains(filename::String)

Loads a protein structure (represented as a `Vector{Protein.Chain}`) from a PDB file.
Assumes that each residue has N, CA, and C atoms.
"""
readchains(filename::AbstractString, format::Type{<:ProteinFileFormat}) = collectchains(read(filename, format))

readpdb(filename::String) = readchains(filename, PDBFormat)
readcif(filename::String) = readchains(filename, MMCIFFormat)

"""
    writepdb(protein::Vector{Protein.Chain}, filename)

Write a protein structure (represented as a `Vector{Protein.Chain}`s) to a PDB file.
"""
function writepdb(protein::Vector{Chain}, filename)
    atom_records = BioStructures.AtomRecord[]
    index = 0
    residue_index = 0
    for chain in protein
        residue_backbone_coords = reshape(chain.backbone.coords, 3, 3, :)
        assign_missing_oxygens!(chain)
        for (i, residue) in enumerate(chain)
            resname = get(threeletter_aa_names, residue.aa, "XXX")
            residue_index += 1
            residue_atoms = [
                [Atom(name, coords) for (name, coords) in zip(BACKBONE_ATOM_NAMES, eachcol(view(residue_backbone_coords, :, :, i)))];
                filter(a -> !(a.name in BACKBONE_ATOM_NAMES), residue.atoms)
            ]
            for atom in residue_atoms
                index += 1
                push!(atom_records, BioStructures.AtomRecord(
                    false,
                    index,
                    atom.name,
                    ' ',
                    resname,
                    chain.id,
                    residue.num,
                    residue.ins_code,
                    atom.coords,
                    1.0,
                    0.0,
                    strip(atom.name)[1:1],
                    "",
                ))
            end
        end
    end
    pdblines = BioStructures.pdbline.(atom_records)
    io = open(filename, "w")
    for line in pdblines
        println(io, line)
    end
    close(io)
    return nothing
end
