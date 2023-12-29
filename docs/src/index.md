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

Backboner is a Julia package that offers a suite of tools for storing protein backbone atom positions, estimating oxygen atom positions, assigning secondary structure, and more. [View the source code on GitHub](https://github.com/MurrellGroup/Backboner.jl) (licensed under MIT).

## Installation

Backboner is a registered Julia package, and can be installed with the Julia package manager:

```julia
using Pkg
Pkg.add("Backboner")
```

## Usage

The `Protein` type wraps a vector of `Chain`s.

```jldoctest
julia> using Backboner, Backboner.Protein

julia> protein = readpdb("test/data/1ZAK.pdb")
2-element Vector{Chain}:
 Chain A with 220 residues
 Chain B with 220 residues

julia> chain = protein["A"] # chains can be accessed by name
Chain A with 220 residues

julia> protein["A"] == protein[1] # numeric indexing also works
true

julia> new_protein = [protein["A"]] # create a new protein with a single chain
1-element Vector{Chain}:
 Chain A with 220 residues

julia> writepdb(new_protein, "test/data/1ZAK_A.pdb");
```