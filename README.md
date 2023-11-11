# Backboner

[![Build Status](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/Backboner.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/Backboner.jl)

Backboner is a Julia package for storing the atom positions of protein backbones in a compact format. The `Protein` type is a wrapper for a vector of `Chain`s, which wraps the `Backbone{4}` type, since it stores the positions of 4 atoms per residue. The `Backbone{N}` type has the `N` type parameter in order to remain flexible. One can load in a protein backbone without the oxygen positions, in the form of a `Backbone{3}` type, and then add the oxygen positions using the `add_oxygen_slice` function.

The package includes a SecondaryStructure type, which describes the secondary structure of a single residue in a chain. The secondary structure of an entire chain is described by a vector of SecondaryStructures. For assignment of secondary structure, this package uses the [AssigningSecondaryStructures.jl](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl) package, which implements a simplified version of the DSSP algorithm to assign three types of secondary structure to residues: loop, helix, and strand.

Protein backbones can be loaded from a PDB file using the `pdb_to_protein` function, which returns a `Protein` object. Inversely, a `Protein` object can be written to a PDB file using the `protein_to_pdb` function.

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