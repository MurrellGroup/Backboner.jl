# Oxygen atoms

The reason for having `N` as a type parameter is to make `Backbone` more flexible by allowing one to pass only the N, CA, and C atoms of a backbone, such that the O atom positions can added in using the `add_oxygens` function.

```jldoctest
julia> using Backboner

julia> protein = pdb_to_protein("test/data/1ZAK.pdb")

julia> chain = protein["A"]

julia> backbone4 = chain.backbone
3×4×220 Backbone{4, Float32}:
[:, :, 1] =
 22.346  22.901  23.227  22.689
 17.547  18.031  16.793  15.72
 23.294  21.993  21.163  21.448

;;; … 

[:, :, 220] =
 21.808  22.263  21.085  19.939
 13.861  13.862  14.233  13.851
  2.734   1.355   0.446   0.791

julia> backbone3 = remove_column(backbone4, 4) # remove oxygen column
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

julia> backbone4_approx = add_oxygens(backbone3) # add oxygen column
3×4×220 Backbone{4, Float32}:
[:, :, 1] =
 22.346  22.901  23.227  22.6697 # decent approximation
 17.547  18.031  16.793  15.7257
 23.294  21.993  21.163  21.4295

;;; … 

[:, :, 220] =
 21.808  22.263  21.085  20.2198 # last oxygen position won't match original
 13.861  13.862  14.233  13.42
  2.734   1.355   0.446   0.112399
```

!!! note
    The `add_oxygens` function adds oxygen atoms to the backbone using idealized geometry, and oxygens atom will on average deviate [0.05 Å](https://github.com/MurrellGroup/Backboner.jl/blob/main/test/backbone/oxygen.jl) from the original positions.
    Moreover, the last oxygen atom is also given a random orientation, as that information is lost when the backbone is reduced to 3 atoms.