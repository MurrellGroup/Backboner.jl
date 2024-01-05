export Residue

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
    aa3 = get(THREE_LETTER_AA_CODES, residue.aa, "XXX")
    ss = residue.ss
    print(io, "$(summary(residue)) $num $aa3 $ss")
end