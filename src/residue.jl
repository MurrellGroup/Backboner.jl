export Residue

struct Residue
    index::Integer
    backbone::Backbone{4}
    aminoacid::Char
    secondarystructure::SecondaryStructure

    function Residue(
        index::Integer,
        backbone::Backbone{4},
        aminoacid::Char = 'G',
        secondarystructure::SecondaryStructure = Unassigned
    )
        return new(index, backbone, aminoacid, secondarystructure)
    end
end

function Base.summary(residue::Residue)
    index = lpad(string(residue.index), length(string(length(residue.backbone))))
    code = get(THREE_LETTER_AA_CODES, residue.aminoacid, "XXX")
    ss = lpad(string(residue.secondarystructure), 6)
    return "Residue $index $code $ss"
end
Base.show(io::IO, residue::Residue) = print(io, summary(residue))