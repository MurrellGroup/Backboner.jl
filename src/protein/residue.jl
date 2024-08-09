const threeletter_aa_names = Dict{Char, String}([Char(v) => k for (k, v) in BioStructures.threeletter_to_aa])

# not user-facing

struct Residue
    num::Integer
    atoms::Vector{Atom}
    aa::Char
    ss::Char
    ins_code::Char

    function Residue(num::Integer, atoms::Vector{Atom}, aa::Char='G', ss::Char=' ', ins_code::Char=' ')
        return new(num, atoms, aa, ss, ins_code)
    end
end

function Base.show(io::IO, residue::Residue)
    num = residue.num
    aa = get(threeletter_aa_names, residue.aa, "XXX")
    ss = residue.ss == ' ' ? "" : "($(residue.ss))"
    ins_code = residue.ins_code == ' ' ? "" : residue.ins_code
    print(io, "Residue $num$ins_code $aa with $(length(residue.atoms)) atom$(length(residue.atoms) == 1 ? "" : "s") $ss")
end