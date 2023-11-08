import PDBTools

export load_pdb_backbone

function collect_ncaco(atoms::Vector{PDBTools.Atom})
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

function Chain(atoms::Vector{PDBTools.Atom})
    id = PDBTools.chain(atoms[1])
    @assert length(Set(PDBTools.chain.(atoms))) == 1 "atoms must be from the same chain"
    residues = collect_ncaco(atoms)
    coords = zeros(Float32, (3, 4, length(residues)))
    for (i, residue) in enumerate(residues)
        for (j, atom) in enumerate(residue)
            coords[:, j, i] = [atom.x, atom.y, atom.z]
        end
    end
    chain = Chain(id, coords)
    return chain
end

function Backbone(atoms::Vector{PDBTools.Atom})
    filter!(a -> a.name in ["N", "CA", "C", "O"], atoms)
    ids = PDBTools.chain.(atoms)
    chains = [Chain(atoms[ids .== id]) for id in unique(ids)]
    backbone = Backbone(chains)
    return backbone
end

"""
    load_pdb_backbones(filename::String)

Assumes that each residue starts with four atoms: N, CA, C, O.
"""
load_pdb_backbone(filename::String) = Backbone(PDBTools.readPDB(filename))