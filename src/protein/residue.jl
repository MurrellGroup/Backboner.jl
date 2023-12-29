export Residue

struct Residue
    index::Integer
    aa::Char
    ss::Char

    function Residue(index::Integer, aa::Char = 'G', ss::Char = ' ')
        return new(index, aa, ss)
    end
end

function Base.summary(residue::Residue)
    index = residue.index
    aa3 = get(THREE_LETTER_AA_CODES, residue.aa, "XXX")
    ss = residue.ss
    return "Residue $index $aa3 $ss"
end

Base.show(io::IO, residue::Residue) = print(io, summary(residue))