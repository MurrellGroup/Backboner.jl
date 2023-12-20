export Residue

struct Residue
    index::Integer
    backbone::Backbone
    aa::Char
    ss::Char

    function Residue(
        index::Integer,
        backbone::Backbone,
        aa::Char = 'G',
        ss::Char = ' ',
    )
        return new(index, backbone, aa, ss)
    end
end

function Base.summary(residue::Residue)
    index = lpad(string(residue.index), ndigits(length(residue.backbone) ÷ 3))
    aa3 = get(THREE_LETTER_AA_CODES, residue.aa, "XXX")
    ss = residue.ss
    return "Residue $index $aa3 $ss"
end

Base.show(io::IO, residue::Residue) = print(io, summary(residue))