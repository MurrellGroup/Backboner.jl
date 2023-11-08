# Backboner

[![Build Status](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/Backboner.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/Backboner.jl)

Backboner is a Julia package for storing the atom positions of protein backbones in a compact format. The `Backbone` type is a wrapper for a vector of `Chain`s. Both of these types have a type parameter `A` for specifying the number of atoms per residue, which is linked to the size of a chain's coordinate array: (3, A, len).

The package includes a SecondaryStructure type, which describes the secondary structure of a single residue in a chain. The secondary structure of an entire chain is thus described by a vector of SecondaryStructures. This package uses the [AssigningSecondaryStructures.jl](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl) package, which implements a simplified version of the DSSP algorithm to assign three types of secondary structure to residues: loop, helix, and strand.

Protein backbones can be loaded from a PDB file using the `load_pdb_backbone` function, which returns a `Backbone` object.