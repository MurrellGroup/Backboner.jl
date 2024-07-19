export readchains, writechains

export readpdb, writepdb, readcif, writecif
#@deprecate readpdb(path::AbstractString) readchains(path::AbstractString)
#@deprecate writepdb(path::AbstractString, chains) writechains(path::AbstractString, chains)

using BioStructures: PDBFormat, MMCIFFormat
const ProteinFileFormat = Union{PDBFormat, MMCIFFormat}
export PDBFormat, MMCIFFormat

#= TODO:
Vector{Chain} to MolecularStructure conversion to allow for writing to other formats than PDB
temporary solution could be to write pdb, read+delete pdb, write e.g. cif
=#

import BioStructures

const AMINOACIDS = Set("ACDEFGHIKLMNPQRSTVWY")
const BACKBONE_ATOM_NAMES = ["N", "CA", "C"]

threeletter_to_oneletter(threeletter::AbstractString) = Char(get(BioStructures.threeletter_to_aa, threeletter, 'X'))
oneletter_resname(residue::BioStructures.AbstractResidue) = threeletter_to_oneletter(BioStructures.resname(residue))

backbone_atom_selector(atom::BioStructures.AbstractAtom) = BioStructures.atomnameselector(atom, BACKBONE_ATOM_NAMES)
backbone_residue_selector(residue::BioStructures.AbstractResidue) =
    oneletter_resname(residue) in AMINOACIDS &&
    BioStructures.countatoms(residue, backbone_atom_selector) == 3 &&
    BioStructures.standardselector(residue) &&
    !BioStructures.disorderselector(residue)

function get_atom(residue::BioStructures.AbstractResidue, atom_name::AbstractString)
    selector = atom -> BioStructures.atomnameselector(atom, [atom_name])
    residue_atoms = BioStructures.collectatoms(residue, selector)
    return only(residue_atoms)
end

function Backboner.Backbone(residue::BioStructures.AbstractResidue)
    atom_name_to_atom_coords = atom_name -> BioStructures.coords(get_atom(residue, atom_name))
    residue_coords = stack(atom_name_to_atom_coords, BACKBONE_ATOM_NAMES)
    return Backbone(residue_coords)
end

function Backboner.Backbone(residues::Vector{BioStructures.AbstractResidue})
    chain_coords = mapreduce(residue -> Backboner.Backbone(residue).coords, hcat, residues; init=Matrix{Float64}(undef, 3, 0))
    return Backbone(chain_coords)
end

Backboner.Backbone(chain::BioStructures.Chain, selector) = Backboner.Backbone(BioStructures.collectresidues(chain, selector))

aminoacid_sequence(residues::Vector{BioStructures.AbstractResidue}) = oneletter_resname.(residues)
aminoacid_sequence(chain::BioStructures.Chain, selector) = aminoacid_sequence(BioStructures.collectresidues(chain, selector))

function get_residue_atoms(residues::Vector{BioStructures.AbstractResidue})
    residue_atoms = Vector{Atom}[]
    for res in residues
        push!(residue_atoms, Atom.(BioStructures.collectatoms(res, a -> BioStructures.standardselector(a) && !BioStructures.disorderselector(a))))
    end
    return residue_atoms
end

function Protein.Chain(residues::Vector{BioStructures.AbstractResidue})
    id = only(unique(BioStructures.chainid.(residues)))
    backbone = Backbone(residues)
    aavector = aminoacid_sequence(residues)
    resnums = BioStructures.resnumber.(residues)
    ins_codes = BioStructures.inscode.(residues)
    modelnum = only(unique(BioStructures.modelnumber.(residues)))
    residue_atoms = get_residue_atoms(residues)
    return Protein.Chain(backbone; id, resnums, ins_codes, aavector, modelnum, residue_atoms)
end

function Protein.Chain(chain::BioStructures.Chain, selector=backbone_residue_selector)
    residues = BioStructures.collectresidues(chain, selector)
    isempty(residues) && return Protein.Chain(BioStructures.chainid(chain), Backbone(Matrix{Float64}(undef, 3, 0)); modelnum=BioStructures.modelnumber(chain))
    return Protein.Chain(residues)
end

function collectchains(struc::BioStructures.MolecularStructure, selector=backbone_residue_selector)
    chains = Protein.Chain[]
    for model in struc, chain in model
        isempty(chain) || push!(chains, Protein.Chain(chain, selector))
    end
    return chains
end

const pdbextension_to_format = Dict(ext => format for (format, ext) in BioStructures.pdbextension)

get_format(path::AbstractString) = get(pdbextension_to_format, lowercase(last(splitext(path))[2:end]), PDBFormat)

"""
    readchains(path, format) -> chains::Vector{Protein.Chain}

Loads a protein structure from a PDB file.

Exported formats: `PDBFormat`, `MMCIFFormat`

## Examples

```jldoctest
julia> readchains("example.pdb") # detects PDB format from extension

julia> readchains("example.cif") # detects mmCIF format from extension

julia> readchains("example.abc", PDBFormat) # force PDB format

julia> readchains("example.xyz", MMCIFFormat) # force mmCIF format
```
"""
readchains(path::AbstractString, format::Type{<:ProteinFileFormat}) = collectchains(read(path, format))
readchains(path::AbstractString) = readchains(path, get_format(path))

"""
    readpdb(path) -> chains::Vector{Protein.Chain}
"""
readpdb(path::AbstractString) = readchains(path, PDBFormat)

"""
    readcif(path) -> chains::Vector{Protein.Chain}
"""
readcif(path::AbstractString) = readchains(path, MMCIFFormat)

"""
    writepdb(path, chains::AbstractVector{Protein.Chain})
    writepdb(path, chain::Protein.Chain)
"""
function writepdb(path::AbstractString, chains::AbstractVector{Protein.Chain})
    atom_records = BioStructures.AtomRecord[]
    index = 0
    residue_index = 0
    for chain in chains
        residue_backbone_coords = reshape(chain.backbone.coords, 3, 3, :)
        assign_missing_oxygens!(chain)
        for (i, residue) in enumerate(chain)
            resname = get(threeletter_aa_names, residue.aa, "XXX") # threletter_aa_names in residue.jl
            residue_index += 1
            residue_atoms = [
                [Atom(name, coords) for (name, coords) in zip(BACKBONE_ATOM_NAMES, eachcol(view(residue_backbone_coords, :, :, i)))];
                filter(a -> !(a.name in BACKBONE_ATOM_NAMES), residue.atoms)]
            for atom in residue_atoms
                index += 1
                push!(atom_records, BioStructures.AtomRecord(false, index, atom.name, ' ', resname, chain.id,
                    residue.num, residue.ins_code, atom.coords, 1.0, 0.0, strip(atom.name)[1:1], ""))
            end
        end
    end
    pdblines = BioStructures.pdbline.(atom_records)
    open(path, "w") do io
        for line in pdblines
            println(io, line)
        end
    end
    return nothing
end

writepdb(path::AbstractString, chain::Protein.Chain) = writepdb(path, [chain])

# compat
writepdb(chains::Union{Protein.Chain, AbstractVector{Protein.Chain}}, path::AbstractString) = writepdb(path, chains)

"""
    writechains(path, chains::AbstractVector{Protein.Chain}, format)
    writechains(path, chain::Protein.Chain, format)

Write a protein structure (represented as a `Vector{Protein.Chain}`s) to file with the specified format.

Exported formats: `PDBFormat`, `MMCIFFormat`

## Examples

```jldoctest
julia> writechains("example.pdb", chains) # detects PDB format from extension

julia> writechains("example.cif", chains) # detects mmCIF format from extension

julia> writechains("example.abc", chains, PDBFormat) # force PDB format

julia> writechains("example.xyz", chains, MMCIFFormat) # force mmCIF format
```
"""
writechains(path::AbstractString, chains::AbstractVector{Protein.Chain}, ::Type{PDBFormat}) = writepdb(path, chains)

function writechains(path::AbstractString, chains::AbstractVector{Protein.Chain}, format::Type{<:ProteinFileFormat})
    struc = mktempdir() do temp_dir
        temp_path = joinpath(temp_dir, "temp.pdb")
        writechains(temp_path, chains, PDBFormat)
        read(temp_path, PDBFormat) # loads BioStructures.MolecularStructure
    end
    write_function = format == PDBFormat ? BioStructures.writepdb : BioStructures.writemmcif
    write_function(path, struc)
end

writechains(path::AbstractString, chains::AbstractVector{Protein.Chain}) = writechains(path, chains, get_format(path))

writechains(path, chain::Protein.Chain, args...) = writechains(path, [chain], args...)