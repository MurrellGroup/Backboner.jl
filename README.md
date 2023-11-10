# Backboner

[![Build Status](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/Backboner.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/Backboner.jl)

Backboner is a Julia package for storing the atom positions of protein backbones in a compact format. The `Protein` type is a wrapper for a vector of `Chain`s, which wraps the `Backbone{4}` type, since it stores the positions of 4 atoms per residue. The `Backbone{N}` type has the `N` type parameter in order to remain flexible. One can load in a protein backbone without the oxygen positions, in the form of a `Backbone{3}` type, and then add the oxygen positions using the `add_oxygen_slice` function.

The package includes a SecondaryStructure type, which describes the secondary structure of a single residue in a chain. The secondary structure of an entire chain is described by a vector of SecondaryStructures. For assignment of secondary structure, this package uses the [AssigningSecondaryStructures.jl](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl) package, which implements a simplified version of the DSSP algorithm to assign three types of secondary structure to residues: loop, helix, and strand.

Protein backbones can be loaded from a PDB file using the `load_pdb` function, which returns a `Protein` object.