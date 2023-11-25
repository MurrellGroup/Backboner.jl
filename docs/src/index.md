```@meta
CurrentModule = VectorizedKmers
DocTestSetup = quote
    using VectorizedKmers
end
```

# Backboner

[![Latest Release](https://img.shields.io/github/release/MurrellGroup/Backboner.jl.svg)](https://github.com/MurrellGroup/Backboner.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://MurrellGroup.github.io/Backboner.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://MurrellGroup.github.io/Backboner.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/Backboner.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/Backboner.jl)

Backboner is a Julia package for storing the atom positions of protein backbones in a compact format, with additional functionality for estimating oxygen atom positions, assigning secondary structure, and loading backbones from residue locations and rotation matrices.

## Types & functions

The `Protein` type wraps a vector of `Chain`s, which in turn wraps the `Backbone{4}` type (4, because it stores the positions of 4 atoms per residue: N, CA, C, O). The `Backbone{N}` type has the `N` type parameter in order to remain flexible. It allows one pass only the N, CA, and C atoms of a backbone, such that the O atom positions can added in using the `add_oxygens` function.

The package includes a `SecondaryStructure` type, which describes the secondary structure of a single residue in a chain. The secondary structure of an entire chain is described by a `Vector{SecondaryStructure}`. For assignment of secondary structure, this package uses the [AssigningSecondaryStructure.jl](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl) package, which implements a simplified version of the DSSP algorithm to assign three types of secondary structure to residues: loop, helix, and strand.

Protein backbones can be loaded from a PDB file using the `pdb_to_protein` function, which returns a `Protein` instance. Inversely, a `Protein` instance can be written to a PDB file using the `protein_to_pdb` function.

## Example

```julia
julia> using Backboner

julia> protein = pdb_to_protein("test/data/1ZAK.pdb")
2-element Protein{Float32}:
 Chain A with 220 residues
 Chain B with 220 residues

julia> chain_A = protein[1]
Chain A with 220 residues

julia> chain_A.backbone
3×4×220 Backbone{4, Float32}:
[:, :, 1] =
 22.346  22.901  23.227  22.689
 17.547  18.031  16.793  15.72
 23.294  21.993  21.163  21.448

[:, :, 2] =
 24.115  24.478  25.289  25.266
 16.923  15.779  14.65   13.511
 20.175  19.336  20.009  19.53

;;; … 

[:, :, 219] =
 22.572  21.48   22.041  22.749
 14.235  14.668  14.866  15.845
  5.844   4.974   3.569   3.29

[:, :, 220] =
 21.808  22.263  21.085  19.939
 13.861  13.862  14.233  13.851
  2.734   1.355   0.446   0.791
```