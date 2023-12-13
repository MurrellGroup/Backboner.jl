export pdb_to_protein, protein_to_pdb

const THREE_LETTER_AA_CODES = Dict(
    'A' => "ALA", 'R' => "ARG", 'N' => "ASN", 'D' => "ASP",
    'C' => "CYS", 'Q' => "GLN", 'E' => "GLU", 'G' => "GLY",
    'H' => "HIS", 'I' => "ILE", 'L' => "LEU", 'K' => "LYS",
    'M' => "MET", 'F' => "PHE", 'P' => "PRO", 'S' => "SER",
    'T' => "THR", 'W' => "TRP", 'Y' => "TYR", 'V' => "VAL",
)

const ONE_LETTER_AA_CODES = Dict(v => k for (k, v) in THREE_LETTER_AA_CODES)

function collect_backbone_atoms(atoms::Vector{PDBTools.Atom})
    residues = Vector{PDBTools.Atom}[]
    i = 1
    while i <= length(atoms) - 3 # Ensure there are at least four atoms left to process
        # Check if the next four atoms are N, CA, C, O in order
        if atoms[i].name == "N" && atoms[i+1].name == "CA" && atoms[i+2].name == "C" && atoms[i+3].name == "O" &&
                all(==(PDBTools.resnum(atoms[i])), PDBTools.resnum.(atoms[i+1:i+3]))
            push!(residues, atoms[i:i+3])
            i += 4
        else
            i += 1
        end
    end
    return residues
end

function Backbone(atoms::Vector{PDBTools.Atom})
    backbone_atoms = collect_backbone_atoms(atoms)
    coords = zeros(Float32, (3, 4, length(backbone_atoms)))
    for (i, residue) in enumerate(backbone_atoms)
        for (j, atom) in enumerate(residue)
            coords[:, j, i] = [atom.x, atom.y, atom.z]
        end
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

function Protein(atoms::Vector{PDBTools.Atom})
    filter!(a -> a.name in ["N", "CA", "C", "O"], atoms)
    ids = PDBTools.chain.(atoms)
    chains = [Chain(atoms[ids .== id]) for id in unique(ids)]
    return Protein(chains)
end

"""
    pdb_to_protein(filename::String)

Assumes that each residue starts with four atoms: N, CA, C, O.
"""
pdb_to_protein(filename::String) = Protein(PDBTools.readPDB(filename))

function protein_to_pdb(protein::Protein, filename, header=:auto, footer=:auto)
    atoms = PDBTools.Atom[]
    index = 0
    residue_index = 0
    for chain in protein
        L = length(chain)
        for (resnum, (residue_coords, aa)) in enumerate(zip(eachslice(chain.backbone.coords, dims=3), chain.aavector))
            resname = get(THREE_LETTER_AA_CODES, aa, "XXX")
            residue_index += 1
            for (name, atom_coord_matrix) in zip(["N", "CA", "C", "O"], eachcol(residue_coords))
                index += 1
                atom = PDBTools.Atom(
                    index = index,
                    name = name,
                    resname = resname,
                    chain = chain.id,
                    resnum = resnum,
                    residue = residue_index,
                    x = atom_coord_matrix[1], y = atom_coord_matrix[2], z = atom_coord_matrix[3],
                )
                push!(atoms, atom)
            end
        end

        # reuses the magic vector from src/backbone/oxygen.jl
        # note: the magic vector is meant to be used to calculate the O atom position,
        # but this is basically using it to get the next N position, 
        # so it's a hacky way to get a "random" OXT atom position
        # not actually random -- it has the same orientation as the next-to-last residue
        index += 1
        last_N, last_CA, last_C, last_O = eachcol(chain.backbone.coords[:, :, L])
        rot_matrix = get_rotation_matrix(last_CA, last_C, last_O)
        OXT_pos = rot_matrix' \ magic_vector + last_C
        OXT_atom = PDBTools.Atom(
            index = index,
            name = "OXT",
            resname = get(THREE_LETTER_AA_CODES, chain.aavector[end], "XXX"),
            chain = chain.id,
            resnum = L,
            residue = residue_index,
            x = OXT_pos[1], y = OXT_pos[2], z = OXT_pos[3],
        )
        push!(atoms, OXT_atom)
    end
    PDBTools.writePDB(atoms, filename, header=header, footer=footer)
end