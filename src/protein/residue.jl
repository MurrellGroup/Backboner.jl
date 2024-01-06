export Residue

const threeletter_aa_names = Dict{Char, String}([Char(v) => k for (k, v) in BioStructures.threeletter_to_aa])

const STANDARD_TRIANGLE_ANGSTROM = [
    -1.066  -0.200   1.266;
     0.645  -0.527  -0.118;
     0.000   0.000   0.000;
] #  N       Ca      C

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
    aa = residue.aa
    ss = residue.ss
    print(io, "Residue $num $aa $ss")
end