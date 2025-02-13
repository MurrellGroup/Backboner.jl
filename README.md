# Backboner

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://MurrellGroup.github.io/Backboner.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://MurrellGroup.github.io/Backboner.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/Backboner.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/Backboner.jl)

This package implements types and functions for working with molecular *backbones*, defined here as continuous chains of bonded atoms.[^1]

Backbones can be represented with different types:
- `Backbone`: a type containing a 3xN matrix of coordinates
- `ChainedBonds`: a type that holds vectors of bond lengths, bond angles, and torsion angles
- `Frames`: a collection of rotations and translations (e.g. for representing orientations and locations of protein residues)

Most functions are implemented especially with differentiability in mind, and can be used in combination with automatic differentiation packages like [Zygote.jl](https://github.com/FluxML/Zygote.jl).

## Installation

```
using Pkg
pkg"add Backboner"
```

## Example usage

```julia
julia> using Backboner

julia> backbone = Backbone(rand(Float32, 3, 9))
10-element Backbone{Float32, Matrix{Float32}}:
 [0.48967552, 0.91008425, 0.5339774]
 [0.2951318, 0.38963223, 0.8952989]
 [0.83763623, 0.5279301, 0.3407849]
 [0.88848716, 0.643387, 0.76827604]
 [0.697279, 0.63588345, 0.0889622]
 [0.08590478, 0.6086006, 0.6478121]
 [0.06308746, 0.6704904, 0.55852276]
 [0.46147835, 0.56259614, 0.7884294]
 [0.9694153, 0.052023113, 0.08127427]

julia> ChainedBonds(backbone)
ChainedBonds{Float32, Vector{Float32}} with 8 bond lengths, 7 bond angles, and 6 torsion angles

julia> is_knotted(backbone)
false
```

[^1]: In some contexts, the term *backbone* may be used more loosely, and allow for atoms that are not part of the main continuous chain of atoms. This package does not support storing e.g. oxygen and beta-carbon atoms in the matrix of coordinates, as they are not part of the continuous chain of atoms.
