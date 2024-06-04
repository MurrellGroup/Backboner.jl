```@meta
CurrentModule = Backboner
DocTestSetup = quote
    using Backboner
end
```

# Backboner

[![Latest Release](https://img.shields.io/github/release/MurrellGroup/Backboner.jl.svg)](https://github.com/MurrellGroup/Backboner.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://MurrellGroup.github.io/Backboner.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/Backboner.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/Backboner.jl)

Backboner is a Julia package that offers a set of types and functions for working with molecular *backbones*: defined here as continuous chains of bonded atoms.[^1] The package provides a few different types for representing backbones:
- `Backbone`: a type containing a 3xN matrix of coordinates
- `ChainedBonds`: a type that holds vectors of bond lengths, bond angles, and dihedral angles
- `Frames`: a collection of rotations and translations (e.g. for representing orientations and locations of protein residues)

The `Protein` submodule contains functions and types for working specifically with proteins. A protein can be loaded from a PDB file using the `Backboner.Protein.readpdb` function, which returns a `Vector{Backboner.Protein.Chain}`. Conversely, a `Vector{Backboner.Protein.Chain}` instance can be written to a PDB file using the `writepdb` function.

[View the source code on GitHub](https://github.com/MurrellGroup/Backboner.jl) (licensed under MIT).

## Installation

Backboner is registered, and can be installed in the Julia REPL. Press `]` to enter pkg mode, and then run:

```
add Backboner
```

## Example usage

```julia
julia> using Backboner, Backboner.Protein

julia> chains = readpdb("test/data/1ZAK.pdb")
2-element Vector{Chain}:
 Chain A with 220 residues
 Chain B with 220 residues

julia> backbone = chains[1].backbone
660-element Backbone{Float32, Matrix{Float32}}:
 [22.346, 17.547, 23.294]
 [22.901, 18.031, 21.993]
 [23.227, 16.793, 21.163]
 [24.115, 16.923, 20.175]
 [24.478, 15.779, 19.336]
 ⋮
 [21.48, 14.668, 4.974]
 [22.041, 14.866, 3.569]
 [21.808, 13.861, 2.734]
 [22.263, 13.862, 1.355]
 [21.085, 14.233, 0.446]

julia> ChainedBonds(backbone)
ChainedBonds{Float32, Vector{Float32}} with 659 bonds, 658 angles, and 657 dihedrals

julia> is_knotted(backbone)
false

julia> import Zygote # unlock the `idealize` method for backbones

julia> idealize(backbone, Float32[1.46, 1.52, 1.33], Float32[1.94, 2.04, 2.13])
660-element Backbone{Float32, Matrix{Float32}}:
 [22.348574, 17.582397, 23.289886]
 [22.90583, 17.977451, 21.999538]
 [23.216103, 16.762234, 21.140835]
 [24.204561, 16.88515, 20.259508]
 [24.52946, 15.827013, 19.307465]
 ⋮
 [21.501173, 14.705252, 4.9825864]
 [22.007494, 14.864742, 3.5582967]
 [21.822643, 13.836198, 2.7356021]
 [22.24875, 13.874594, 1.3396943]
 [21.091076, 14.233609, 0.42247167]
```

[^1]: In some contexts, the term *backbone* may be used more loosely, and allow for atoms that are not part of the main continuous chain of atoms. This package does not support storing e.g. oxygen and beta-carbon atoms in the matrix of coordinates, as they are not part of the continuous chain of atoms.
