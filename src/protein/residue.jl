export Residue

const threeletter_aa_names = Dict{Char, String}([Char(v) => k for (k, v) in BioStructures.threeletter_to_aa])

struct Residue
    num::Integer
    aa::Char
    ss::Char

    function Residue(num::Integer, aa::Char='G', ss::Char=' ')
        return new(num, aa, ss)
    end
end

function Base.show(io::IO, residue::Residue)
    num = residue.num
    aa = get(threeletter_aa_names, residue.aa, "XXX")
    ss = residue.ss
    print(io, "Residue $num $aa $ss")
end