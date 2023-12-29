module Backboner

using LinearAlgebra

include("backbone.jl")
include("bonds.jl")
include("protein/protein.jl")

using .Protein

pdb_to_protein, protein_to_pdb = readpdb, writepdb

export pdb_to_protein, protein_to_pdb, readpdb, writepdb


end