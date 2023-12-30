export readpdb, writepdb, pdb_to_protein, protein_to_pdb

const THREE_LETTER_AA_CODES = Dict(
    'A' => "ALA", 'R' => "ARG", 'N' => "ASN", 'D' => "ASP",
    'C' => "CYS", 'Q' => "GLN", 'E' => "GLU", 'G' => "GLY",
    'H' => "HIS", 'I' => "ILE", 'L' => "LEU", 'K' => "LYS",
    'M' => "MET", 'F' => "PHE", 'P' => "PRO", 'S' => "SER",
    'T' => "THR", 'W' => "TRP", 'Y' => "TYR", 'V' => "VAL",
)

const ONE_LETTER_AA_CODES = Dict(v => k for (k, v) in THREE_LETTER_AA_CODES)

function collect_backbone_atoms(atoms::Vector{PDBTools.Atom})
    backbone_atoms = PDBTools.Atom[]
    i = 1
    while i <= length(atoms) - 2 # Ensure there are at least three atoms left to process
        # Check if the next four atoms are N, CA, C, O in order
        if atoms[i].name == "N" && atoms[i+1].name == "CA" && atoms[i+2].name == "C" && allequal(PDBTools.resnum.(atoms[i:i+2]))
            append!(backbone_atoms, atoms[i:i+2])
            i += 3
        else
            i += 1
        end
    end
    return backbone_atoms
end

function Backboner.Backbone(atoms::Vector{PDBTools.Atom})
    backbone_atoms = collect_backbone_atoms(atoms)
    coords = zeros(Float32, (3, length(backbone_atoms)))
    for (i, atom) in enumerate(backbone_atoms)
        coords[:, i] = [atom.x, atom.y, atom.z]
    end
    return Backbone(coords)
end

function Chain(atoms::Vector{PDBTools.Atom})
    id = PDBTools.chain(atoms[1])
    @assert allequal(PDBTools.chain.(atoms)) "atoms must be from the same chain"
    backbone = Backbone(atoms)
    aavector = [get(ONE_LETTER_AA_CODES, atom.resname, 'X') for atom in atoms if atom.name == "CA"]
    return Chain(id, backbone, aavector=aavector)
end

"""
    readpdb(filename::String)

Loads a protein (represented as a `Vector{Protein.Chain}`) from a PDB file.
Assumes that each residue starts with three atoms: N, CA, C.
"""
function readpdb(filename::String)
    atoms = PDBTools.readPDB(filename)
    filter!(a -> a.name in ["N", "CA", "C"], atoms)
    ids = PDBTools.chain.(atoms)
    chains = [Chain(atoms[ids .== id]) for id in unique(ids)]
    return chains
end

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
        for (resnum, (residue_coords, aa)) in enumerate(zip(eachslice(coords, dims=3), chain.aavector))
            resname = get(THREE_LETTER_AA_CODES, aa, "XXX")
            residue_index += 1
            for (name, atom_coords) in zip(["N", "CA", "C", "O"], eachcol(residue_coords))
                index += 1
                atom = PDBTools.Atom(
                    index = index,
                    name = name,
                    resname = resname,
                    chain = chain.id,
                    resnum = resnum,
                    residue = residue_index,
                    x = atom_coords[1],
                    y = atom_coords[2],
                    z = atom_coords[3],
                )
                push!(atoms, atom)
            end
        end
    end
    PDBTools.writePDB(atoms, filename, header=header, footer=footer)
end

pdb_to_protein, protein_to_pdb = readpdb, writepdb