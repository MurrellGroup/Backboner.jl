# Types

## Backbone

The `Backbone` type is designed to efficiently store and manipulate the three-dimensional coordinates of backbone atoms in proteins. It is a fundamental component, and is literally the backbone of protein structures.

`Backbone{N, T}` is a wrapper around a 3xNxL array, where:
- **3** are the three spatial dimensions for the coordinates.
- **N** is the number of atoms per residue.
- **L** is the number of residues in the backbone.
- **T** is the element type of the coordinate array.

## Chain

A `Chain` represents a protein chain, and holds an identifier (usually a single letter), backbone atom coordinates, the amino acid sequence, and secondary structure information.

- `backbone`: An instance of `Backbone{4}`, storing the coordinates of backbone atoms.
- `aavector`: A vector for storing the amino acid sequence.
- `ssvector`: A vector for storing the secondary structure.

The `Chain` type is designed to provide a comprehensive and consistent representation of a protein chain, ensuring that the backbone coordinates align with the corresponding amino acid sequences and secondary structures.

## Protein

The `Protein` type holds multiple `Chain` instances, representing complete protein structures.

- Stores a collection of `Chain` objects.
- Includes a dictionary for quick access to chains via their identifiers.
