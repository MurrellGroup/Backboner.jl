#module Records

import BioStructures

struct GeneralRecord
    id::String
    resnumbers::Vector{Int64}
    chainids::Vector{Int32}
    sequence::String
    rotations::Matrix{Float32}
    locations::Matrix{Float32}
    coordinates::Array{Float32, 3}

    function GeneralRecord(
        id::String,
        resnumbers::Vector{Int64},
        chainids::Vector{Int32},
        sequence::String,
        rotations::Matrix{Float32},
        locations::Matrix{Float32},
        coordinates::Array{Float32, 3},
    )
        L = length(resnumbers)
        length(chainids) == L || throw(ArgumentError("chainids array size mismatch"))
        length(sequence) == L || throw(ArgumentError("sequence array size mismatch"))
        size(rotations, 2) == L || throw(ArgumentError("rotations array size mismatch"))
        size(locations, 2) == L || throw(ArgumentError("locations array size mismatch"))
        size(coordinates, 3) == L || throw(ArgumentError("coordinates array size mismatch"))
        new(id, resnumbers, chainids, sequence, rotations, locations, coordinates)
    end
end

function GeneralRecord(pdbfile::String; dir = ".", id = first(splitext(pdbfile)))
    struc = read(joinpath(dir, pdbfile), BioStructures.PDB)
    seqres = read_chainwise_seqres(joinpath(dir, pdbfile))

    frames = Frames(struc)

    return GeneralRecord(
        id,
        resnumbers(struc, seqres),
        chainids(struc),
        sequence(struc),
        frames.rotations,
        frames.locations ./ 10,
        reshape(Backbone(struc).coords, 3, 3, :) ./ 10
    )
end

# end