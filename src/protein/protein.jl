module Protein

using ..Backboner
import BioStructures

include("residue.jl")
include("chain.jl")
include("oxygen.jl")
include("pdb.jl")
include("idealization.jl")

using PrecompileTools

@compile_workload begin
    chains = readpdb("test/data/1ASS.pdb")
end

end