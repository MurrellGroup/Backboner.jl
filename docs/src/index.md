```@meta
CurrentModule = Backboner
DocTestSetup = quote
    using Backboner
end
```

# Backboner

[![Latest Release](https://img.shields.io/github/release/MurrellGroup/Backboner.jl.svg)](https://github.com/MurrellGroup/Backboner.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://MurrellGroup.github.io/Backboner.jl/stable/)
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
3×660 Backbone{Float32, Matrix{Float32}}:
 22.346  22.901  23.227  24.115  24.478  …  21.48   22.041  21.808  22.263  21.085
 17.547  18.031  16.793  16.923  15.779     14.668  14.866  13.861  13.862  14.233
 23.294  21.993  21.163  20.175  19.336      4.974   3.569   2.734   1.355   0.446

julia> ChainedBonds(backbone)
ChainedBonds{Float32, Vector{Float32}} with 659 bonds, 658 angles, and 657 dihedrals

julia> is_knotted(backbone)
false

julia> import Zygote # unlock the `idealize` method for backbones

julia> idealize(backbone, Float32[1.46, 1.52, 1.33], Float32[1.94, 2.04, 2.13])
3×660 Backbone{Float32, Matrix{Float32}}:
 22.3486  22.9058  23.2161  24.2046  24.5295  …  23.7832   24.2534   23.9791   24.3124   23.1496
 17.5824  17.9775  16.7622  16.8852  15.827      14.3215   14.1375   12.9715   12.7012   13.0422
 23.2899  21.9995  21.1408  20.2595  19.3075      9.99834   8.56466   7.98674   6.59122   5.67358
```

[^1]: In some contexts, the term *backbone* may be used more loosely, and allow for atoms that are not part of the main continuous chain of atoms. This package does not support storing e.g. oxygen and beta-carbon atoms in the matrix of coordinates, as they are not part of the continuous chain of atoms.