export Residue, backbone_atom_coord_matrix

struct Residue
    index::Integer
    backbone::Backbone{4}
    aa::Char
    ss::Char

    function Residue(
        index::Integer,
        backbone::Backbone{4},
        aa::Char = 'G',
        ss::Char = ' ',
    )
        return new(index, backbone, aa, ss)
    end
end

backbone_atom_coord_matrix(residue::Residue) = residue.backbone[residue.index]

function Base.summary(residue::Residue)
    index = lpad(string(residue.index), length(string(length(residue.backbone))))
    aa3 = get(THREE_LETTER_AA_CODES, residue.aa, "XXX")
    ss = residue.ss
    return "Residue $index $aa3 $ss"
end

Base.show(io::IO, residue::Residue) = print(io, summary(residue))