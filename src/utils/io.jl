export pdb_to_protein, protein_to_pdb

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
    return Chain(id, backbone)
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

function protein_to_pdb(protein::Protein, filename::String)
    atoms = PDBTools.Atom[]
    index = 0
    residue_index = 0
    for chain in protein
        backbone = chain.backbone
        for (resnum, residue_coords) in enumerate(eachslice(backbone.coords, dims=3))
            resname = "GLY"
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
    PDBTools.writePDB(atoms, filename)
end