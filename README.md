# Backboner

[![Latest Release](https://img.shields.io/github/release/MurrellGroup/Backboner.jl.svg)](https://github.com/MurrellGroup/Backboner.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://MurrellGroup.github.io/Backboner.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://MurrellGroup.github.io/Backboner.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/Backboner.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/Backboner.jl)

Backboner is a Julia package that offers a set of types and functions for working with molecular *backbones*: defined here as continuous chains of bonded atoms.[^1] The package provides a variety of different types for representing backbones:
- `Backbone`: a type containing a 3xN matrix of coordinates
- `ChainedBonds`: a type that holds vectors of bond lengths, bond angles, and dihedral angles
- `Frames`: a collection of rotations and translations (e.g. for representing orientations and locations of protein residues)

It also includes functions for working with protein chains, and reading/writing PDB files.

## Installation

Backboner is registered, and can be installed with the Julia package manager:

```julia
using Pkg
Pkg.add("Backboner")
```


## Overview

The `Backbone` type is a wrapper for a 3xN matrix of coordinates representing absolute positions of atoms of a continuous molecular backbone. For working with the geometry of a backbone, the `ChainedBonds` type exists to store bond lengths, bond angles, and dihedral angles of a continuous chain of bonds.

The `Protein` submodule contains functions and types for working specifically with proteins. A protein can be loaded from a PDB file using the `Backboner.Protein.readpdb` function, which returns a `Vector{Backboner.Protein.Chain}`. Inversely, a `Vector{Backboner.Protein.Chain}` instance can be written to a PDB file using the `writepdb` function.

## Example

```julia
julia> using Backboner, Backboner.Protein

julia> protein = readpdb("test/data/1ZAK.pdb")
2-element Vector{Chain}:
 Chain A with 220 residues
 Chain B with 220 residues

julia> chain = protein["A"]
Chain A with 220 residues

julia> chain.backbone
660-element Backbone{Float32, Matrix{Float32}}:
 [22.346, 17.547, 23.294]
 [22.901, 18.031, 21.993]
 [23.227, 16.793, 21.163]
 [24.115, 16.923, 20.175]
 ⋮
 [22.041, 14.866, 3.569]
 [21.808, 13.861, 2.734]
 [22.263, 13.862, 1.355]
 [21.085, 14.233, 0.446]
```

[^1]: In some contexts, the term *backbone* may be used more loosely, and allow for atoms that do not part of the main continuous chain of atoms. This package does not support storing e.g. oxygen and beta-carbon atoms in the matrix of coordinates, as they are not part of the continuous chain of atoms.