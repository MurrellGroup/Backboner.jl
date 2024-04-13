var documenterSearchIndex = {"docs":
[{"location":"backboner/#Backboner-API","page":"Backboner API","title":"Backboner API","text":"","category":"section"},{"location":"backboner/","page":"Backboner API","title":"Backboner API","text":"Modules = [Backboner]","category":"page"},{"location":"backboner/#Backboner.Backbone","page":"Backboner API","title":"Backboner.Backbone","text":"Backbone{T <: Real} <: AbstractMatrix{T}\n\nThe Backbone type is designed to efficiently store and manipulate the three-dimensional coordinates of backbone atoms.\n\nExamples\n\nA Backbone can be created from a matrix of coordinates:\n\njulia> backbone = Backbone(zeros(3, 5)) # 5 atoms with 3 coordinates each\n3×5 Backbone{Float64, Matrix{Float64}}:\n 0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0\n\njulia> backbone[1] = [1.0, 2.0, 3.0]; # set the first atom's coordinates\n\njulia> backbone\n3×5 Backbone{Float64, Matrix{Float64}}:\n 1.0  0.0  0.0  0.0  0.0\n 2.0  0.0  0.0  0.0  0.0\n 3.0  0.0  0.0  0.0  0.0\n\njulia> backbone[1:2] # indexing by range returns a new Backbone\n3×2 Backbone{Float64, Matrix{Float64}}:\n 1.0  0.0\n 2.0  0.0\n 3.0  0.0\n\nArrays will always be flattened to a 3xN matrix:\n\njulia> backbone = Backbone(zeros(3, 3, 4)) # e.g. 3 coordinates per atom, 3 atoms per residue, 4 residues\n3×12 Backbone{Float64, Matrix{Float64}}:\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n","category":"type"},{"location":"backboner/#Backboner.ChainedBonds","page":"Backboner API","title":"Backboner.ChainedBonds","text":"ChainedBonds{T <: Real, V <: AbstractVector{T}}\n\nA lazy way to store a backbone as a series of bond lengths, angles, and dihedrals.\n\nExamples\n\njulia> backbone = Protein.readpdb(\"test/data/1ZAK.pdb\")[\"A\"].backbone\n3×660 Backbone{Float32, Matrix{Float32}}:\n 22.346  22.901  23.227  24.115  24.478  25.289  26.091  26.814  …  23.137  22.572  21.48   22.041  21.808  22.263  21.085\n 17.547  18.031  16.793  16.923  15.779  14.65   14.958  13.827     13.041  14.235  14.668  14.866  13.861  13.862  14.233\n 23.294  21.993  21.163  20.175  19.336  20.009  21.056  21.652      5.676   5.844   4.974   3.569   2.734   1.355   0.446\n\njulia> bonds = ChainedBonds(backbone)\nChainedBonds{Float32, Vector{Float32}} with 659 bonds, 658 angles, and 657 dihedrals\n\n\n\n\n\n","category":"type"},{"location":"backboner/#Backboner.Frame","page":"Backboner API","title":"Backboner.Frame","text":"Frame{T <: Real}\n\nA Frame is a combination of a rotation and a translation, which can be applied to a set of coordinates.\n\n\n\n\n\n","category":"type"},{"location":"backboner/#Backboner.Frames","page":"Backboner API","title":"Backboner.Frames","text":"Frames{T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{Frame{T}}\n\nThe Frames type is designed to efficiently store and manipulate the rotation and translation of a set of Frames.\n\n\n\n\n\n","category":"type"},{"location":"backboner/#Backboner.append_bonds-Tuple{Backbone, AbstractVector{<:Real}, AbstractVector{<:Real}, AbstractVector{<:Real}}","page":"Backboner API","title":"Backboner.append_bonds","text":"append_bonds(backbone, lengths, angles, dihedrals)\n\n\n\n\n\n","category":"method"},{"location":"backboner/#Backboner.idealize","page":"Backboner API","title":"Backboner.idealize","text":"note: Note\nZygote must be imported in order to activate the ZygoteIdealizationExt extension, which defines the idealize(::Backbone) method.\n\n\n\n\n\n","category":"function"},{"location":"backboner/#Backboner.is_knotted","page":"Backboner API","title":"Backboner.is_knotted","text":"is_knotted(backbone::Backbone)\n\nCheck if a backbone is knotted.\n\n\n\n\n\n","category":"function"},{"location":"backboner/#Backboner.prepend_bonds-Tuple{Backbone, AbstractVector{<:Real}, AbstractVector{<:Real}, AbstractVector{<:Real}}","page":"Backboner API","title":"Backboner.prepend_bonds","text":"prepend_bonds(backbone, lengths, angles, dihedrals)\n\n\n\n\n\n","category":"method"},{"location":"","page":"Overview","title":"Overview","text":"CurrentModule = Backboner\nDocTestSetup = quote\n    using Backboner\nend","category":"page"},{"location":"#Backboner","page":"Overview","title":"Backboner","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"(Image: Latest Release) (Image: MIT license) (Image: Documentation) (Image: Documentation) (Image: Build Status) (Image: Coverage)","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"Backboner is a Julia package that offers a set of types and functions for working with molecular backbones: defined here as continuous chains of bonded atoms.[1] The package provides a few different types for representing backbones:","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"Backbone: a type containing a 3xN matrix of coordinates\nChainedBonds: a type that holds vectors of bond lengths, bond angles, and dihedral angles\nFrames: a collection of rotations and translations (e.g. for representing orientations and locations of protein residues)","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"The Protein submodule contains functions and types for working specifically with proteins. A protein can be loaded from a PDB file using the Backboner.Protein.readpdb function, which returns a Vector{Backboner.Protein.Chain}. Conversely, a Vector{Backboner.Protein.Chain} instance can be written to a PDB file using the writepdb function.","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"View the source code on GitHub (licensed under MIT).","category":"page"},{"location":"#Installation","page":"Overview","title":"Installation","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"Backboner is registered, and can be installed in the Julia REPL. Press ] to enter pkg mode, and then run:","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"add Backboner","category":"page"},{"location":"#Example-usage","page":"Overview","title":"Example usage","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"julia> using Backboner, Backboner.Protein\n\njulia> chains = readpdb(\"test/data/1ZAK.pdb\")\n2-element Vector{Chain}:\n Chain A with 220 residues\n Chain B with 220 residues\n\njulia> backbone = chains[1].backbone\n3×660 Backbone{Float32, Matrix{Float32}}:\n 22.346  22.901  23.227  24.115  24.478  …  21.48   22.041  21.808  22.263  21.085\n 17.547  18.031  16.793  16.923  15.779     14.668  14.866  13.861  13.862  14.233\n 23.294  21.993  21.163  20.175  19.336      4.974   3.569   2.734   1.355   0.446\n\njulia> ChainedBonds(backbone)\nChainedBonds{Float32, Vector{Float32}} with 659 bonds, 658 angles, and 657 dihedrals\n\njulia> is_knotted(backbone)\nfalse\n\njulia> import Zygote # unlock the `idealize` method for backbones\n\njulia> idealize(backbone, Float32[1.46, 1.52, 1.33], Float32[1.94, 2.04, 2.13])\n3×660 Backbone{Float32, Matrix{Float32}}:\n 22.3486  22.9058  23.2161  24.2046  24.5295  …  23.7832   24.2534   23.9791   24.3124   23.1496\n 17.5824  17.9775  16.7622  16.8852  15.827      14.3215   14.1375   12.9715   12.7012   13.0422\n 23.2899  21.9995  21.1408  20.2595  19.3075      9.99834   8.56466   7.98674   6.59122   5.67358","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"[1]: In some contexts, the term backbone may be used more loosely, and allow for atoms that are not part of the main continuous chain of atoms. This package does not support storing e.g. oxygen and beta-carbon atoms in the matrix of coordinates, as they are not part of the continuous chain of atoms.","category":"page"},{"location":"protein/#Protein-API","page":"Protein API","title":"Protein API","text":"","category":"section"},{"location":"protein/","page":"Protein API","title":"Protein API","text":"These following functions need to be imported explicitly with using Backboner.Protein","category":"page"},{"location":"protein/","page":"Protein API","title":"Protein API","text":"Modules = [Backboner.Protein]","category":"page"},{"location":"protein/#Backboner.Protein.Chain","page":"Protein API","title":"Backboner.Protein.Chain","text":"Chain <: AbstractVector{Residue}\n\nA Chain represents a chain of a protein, and is a vector of Residues, which are instantiated from indexing the chain.\n\nFields\n\nid::String: A string identifier (usually a single letter).\nbackbone::Backbone: A backbone with a length divisible by 3, to ensure 3 atoms per residue (N, Ca, C).\nmodelnum::Int: The model number of the chain.\nresnums::Vector{Int}: storing the residue numbers.\naavector::Vector{Char}: storing the amino acid sequence.\nssvector::Vector{Char}: storing the secondary structure.\n\n\n\n\n\n","category":"type"},{"location":"protein/#Backboner.Protein.alphacarbon_angles-Tuple{Backboner.Protein.Chain}","page":"Protein API","title":"Backboner.Protein.alphacarbon_angles","text":"alphacarbon_angles(chain::Chain)\nalphacarbon_angles(bonds::ChainedBonds)\n\nCalculate the angles at the alphacarbon atoms (N-Ca-C angles) of a chains backbone, or take directly from a precalculated ChainedBonds instance.\n\n\n\n\n\n","category":"method"},{"location":"protein/#Backboner.Protein.alphacarbon_carbonyl_distances-Tuple{Backboner.Protein.Chain}","page":"Protein API","title":"Backboner.Protein.alphacarbon_carbonyl_distances","text":"alphacarbon_carbonyl_distances(chain::Chain)\nalphacarbon_carbonyl_distances(backbone::Backbone)\nalphacarbon_carbonyl_distances(bonds::ChainedBonds)\n\nCalculate the distances between all pairs of contiguous alpha-carbon and carbonyl atoms in a chain. Returns a vector of distances of length length(chain).\n\n\n\n\n\n","category":"method"},{"location":"protein/#Backboner.Protein.alphacarbon_coords-Tuple{Backboner.Protein.Chain}","page":"Protein API","title":"Backboner.Protein.alphacarbon_coords","text":"alphacarbon_coords(chain::Chain)\nalphacarbon_coords(backbone::Backbone)\n\nReturns the coordinates of all alphacarbon atoms in a chain, as a 3xN matrix.\n\n\n\n\n\n","category":"method"},{"location":"protein/#Backboner.Protein.carbonyl_angles-Tuple{Backboner.Protein.Chain}","page":"Protein API","title":"Backboner.Protein.carbonyl_angles","text":"carbonyl_angles(chain::Chain)\ncarbonyl_angles(bonds::ChainedBonds)\n\nCalculate the angles at the carbonyl atoms (Ca-C-N angles) of a chain's backbone, or take directly from a precalculated ChainedBonds instance.\n\n\n\n\n\n","category":"method"},{"location":"protein/#Backboner.Protein.carbonyl_coords-Tuple{Backboner.Protein.Chain}","page":"Protein API","title":"Backboner.Protein.carbonyl_coords","text":"carbonyl_coords(chain::Chain)\ncarbonyl_coords(backbone::Backbone)\n\nReturns the coordinates of all carbonyl atoms in a chain, as a 3xN matrix.\n\n\n\n\n\n","category":"method"},{"location":"protein/#Backboner.Protein.carbonyl_nitrogen_distances-Tuple{Backboner.Protein.Chain}","page":"Protein API","title":"Backboner.Protein.carbonyl_nitrogen_distances","text":"carbonyl_nitrogen_distances(chain::Chain)\ncarbonyl_nitrogen_distances(backbone::Backbone)\ncarbonyl_nitrogen_distances(bonds::ChainedBonds)\n\nCalculate the distances between all pairs of contiguous carbonyl and nitrogen atoms in a chain. Returns a vector of distances of length length(chain) - 1.\n\n\n\n\n\n","category":"method"},{"location":"protein/#Backboner.Protein.nitrogen_alphacarbon_distances-Tuple{Backboner.Protein.Chain}","page":"Protein API","title":"Backboner.Protein.nitrogen_alphacarbon_distances","text":"nitrogen_alphacarbon_distances(chain::Chain)\nnitrogen_alphacarbon_distances(backbone::Backbone)\nnitrogen_alphacarbon_distances(bonds::ChainedBonds)\n\nCalculate the distances between all pairs of contiguous nitrogen and alpha-carbon atoms in a chain. Returns a vector of distances of length length(chain).\n\n\n\n\n\n","category":"method"},{"location":"protein/#Backboner.Protein.nitrogen_angles-Tuple{Backboner.Protein.Chain}","page":"Protein API","title":"Backboner.Protein.nitrogen_angles","text":"nitrogen_angles(chain::Chain)\nnitrogen_angles(backbone::Backbone)\nnitrogen_angles(bonds::ChainedBonds)\n\nCalculate the angles at the nitrogen atoms (C-N-Ca angles) of a chains backbone, or take directly from a precalculated ChainedBonds instance.\n\n\n\n\n\n","category":"method"},{"location":"protein/#Backboner.Protein.nitrogen_coords-Tuple{Backboner.Protein.Chain}","page":"Protein API","title":"Backboner.Protein.nitrogen_coords","text":"nitrogen_coords(chain::Chain)\nnitrogen_coords(backbone::Backbone)\n\nReturns the coordinates of all nitrogen atoms in a chain, as a 3xN matrix.\n\n\n\n\n\n","category":"method"},{"location":"protein/#Backboner.Protein.omega_angles-Tuple{Backboner.Protein.Chain}","page":"Protein API","title":"Backboner.Protein.omega_angles","text":"omega_angles(chain::Chain)\nomega_angles(backbone::Backbone)\nomega_angles(bonds::ChainedBonds)\n\nCalculate the omega (Ω) angles of a chain's backbone, or take directly from a precalculated ChainedBonds instance.\n\n\n\n\n\n","category":"method"},{"location":"protein/#Backboner.Protein.oxygen_coords-Tuple{Backboner.Protein.Chain}","page":"Protein API","title":"Backboner.Protein.oxygen_coords","text":"oxygen_coords(chain::Chain)\n\nAdd oxygen atoms to the backbone of a protein, turning the coordinate array from size 3x3xL to 3x4xL-1, where L is the length of the backbone.\n\nExample\n\njulia> chains = readpdb(\"test/data/1ZAK.pdb\")\n2-element Vector{Chain}:\n Chain A with 220 residues\n Chain B with 220 residues\n\njulia> oxygen_coords(chains[\"A\"]) # returns the estimated position of oxygen atoms in chain A (~0.05 Å mean deviation)\n3×220 Matrix{Float32}:\n 22.6697  25.1719  24.7761  25.8559  …  24.7911   22.7649   22.6578   21.24\n 15.7257  13.505   13.5151  11.478      15.0888   12.2361   15.8825   14.2933\n 21.4295  19.5663  22.8638  25.3283      7.95346   4.81901   3.24164  -0.742424        \n\nnote: Note\nThe oxygen_coords function finds the oxygen atoms to the backbone using idealized geometry, and oxygens atom will on average deviate 0.05 Å from original PDB positions. Moreover, the last oxygen atom is essentially given a random (although deterministic) orientation, as that information is lost when the backbone is reduced to 3 atoms, and there's no next nitrogen atom to compare with.\n\n\n\n\n\n","category":"method"},{"location":"protein/#Backboner.Protein.phi_angles-Tuple{Backboner.Protein.Chain}","page":"Protein API","title":"Backboner.Protein.phi_angles","text":"phi_angles(chain::Chain)\nphi_angles(backbone::Backbone)\nphi_angles(bonds::ChainedBonds)\n\nCalculate the phi (φ) angles of a chain's backbone, or take directly from a precalculated ChainedBonds instance.\n\n\n\n\n\n","category":"method"},{"location":"protein/#Backboner.Protein.psi_angles-Tuple{Backboner.Protein.Chain}","page":"Protein API","title":"Backboner.Protein.psi_angles","text":"psi_angles(chain::Chain)\npsi_angles(backbone::Backbone)\npsi_angles(bonds::ChainedBonds)\n\nCalculate the psi (ψ) angles of a chain's backbone, or take directly from a precalculated ChainedBonds instance.\n\n\n\n\n\n","category":"method"},{"location":"protein/#Backboner.Protein.readpdb-Tuple{String}","page":"Protein API","title":"Backboner.Protein.readpdb","text":"readpdb(pdbfile::String)\n\nLoads a protein (represented as a Vector{Protein.Chain}) from a PDB file. Assumes that each residue starts with three atoms: N, CA, C.\n\n\n\n\n\n","category":"method"},{"location":"protein/#Backboner.Protein.writepdb","page":"Protein API","title":"Backboner.Protein.writepdb","text":"writepdb(protein::Vector{Protein.Chain}, filename)\n\nWrite a protein (represented as a Vector{Protein.Chain}s) to a PDB file.\n\n\n\n\n\n","category":"function"}]
}
