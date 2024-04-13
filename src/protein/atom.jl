export Atom

struct Atom
    name::String
    coords::AbstractVector{<:Real}

    function Atom(name::AbstractString, coords::AbstractVector{<:Real})
        length(coords) == 3 || throw(ArgumentError("coords must be a 3-vector"))
        return new(strip(name), coords)
    end
end

Atom(atom::BioStructures.Atom) = Atom(atom.name, atom.coords)

Base.summary(atom::Atom) = "Atom $(atom.name) at $(atom.coords)"
Base.show(io::IO, atom::Atom) = print(io, summary(atom))