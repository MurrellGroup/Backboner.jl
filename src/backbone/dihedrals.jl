export Dihedrals

# source: en.wikipedia.org/wiki/Dihedral_angle#In_polymer_physics
function dihedral_angle(
    u1::T, u2::T, u3::T,
) where T <: AbstractVector
    c12 = cross(u1, u2)
    c23 = cross(u2, u3)
    atan2(
        dot(u2, cross(c12, c23)),
        abs(u2) * dot(c12, c23)
    )
end

struct Dihedrals{T} <: AbstractVector{T}
    angles::AbstractVector{T}

    function Dihedrals(bonds::ContinuousBonds{T}) where T
        angles = Vector{T}(undef, length(bonds)-2)
        vectors = eachcol(bonds.vectors)
        for (i, j) in zip(eachindex(angles), eachindex(vectors))
            angles[i] = dihedral_angle(vectors[j], vectors[j+1], vectors[j+2])
        end
        return new{T}(angles)
    end

    Dihedrals(backbone::Backbone) = Dihedrals(ContinuousBonds(backbone))
end

@inline Base.size(dihedrals::Dihedrals) = size(dihedrals.angles)
@inline Base.length(dihedrals::Dihedrals) = size(dihedrals, 2)
@inline Base.getindex(dihedrals::Dihedrals, i) = dihedrals.angles[i]