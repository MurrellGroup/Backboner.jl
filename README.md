# Backboner

[![Latest Release](https://img.shields.io/github/release/MurrellGroup/Backboner.jl.svg)](https://github.com/MurrellGroup/Backboner.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://MurrellGroup.github.io/Backboner.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://MurrellGroup.github.io/Backboner.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/Backboner.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/Backboner.jl)

Backboner is a Julia package that offers a set of types and functions for working with molecular backbones. It also includes functions for working with protein chains, reading/writing PDB files and assigning secondary structure.

## Installation

Backboner is a registered Julia package, and can be installed with the Julia package manager:

```julia
using Pkg
Pkg.add("Backboner")
```

## Overview

The `Backbone` type is a wrapper for a 3xN matrix of coordinates representing absolute positions of atoms of a continuous molecular backbone. For working with the geometry of a backbone, the `ChainedBonds` type exists to store bond lengths, bond angles, and dihedral angles of a continuous chain of bonds.

The `Protein` submodule contains functions and types for working specifically with proteins. A protein can be loaded from a PDB file using the `Backboner.Protein.readpdb` function, which returns a `Vector{Backboner.Protein.Chain}`. Inversely, a `Vector{Backboner.Protein.Chain}` instance can be written to a PDB file using the `writepdb` function.

**Note:** The `Protein` submodule is not exported by default, but can be imported explicitly with: `using/import Backboner.Protein`.

## Example

```julia
julia> using Backboner

julia> protein = readpdb("test/data/1ZAK.pdb")
2-element Vector{Chain}:
 Chain A with 220 residues
 Chain B with 220 residues

julia> chain = protein["A"]
Chain A with 220 residues

julia> chain.backbone
3×3×220 Backbone{3, Float32}:
[:, :, 1] =
 22.346  22.901  23.227
 17.547  18.031  16.793
 23.294  21.993  21.163

;;; … 

[:, :, 220] =
 21.808  22.263  21.085
 13.861  13.862  14.233
  2.734   1.355   0.446
```
