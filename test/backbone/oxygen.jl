@testset "oxygen.jl" begin

    protein = pdb_to_protein("data/1ASS.pdb")
    backbone4 = protein[1].backbone
    backbone3 = remove_column(backbone4, 4)
    backbone4est = add_oxygens(backbone3)
    true_o_coords = backbone4[1][:, 4, 1:size(backbone4[1], 3)-1] # backbone 4 gets 1 residue shorter
    est_o_coords = backbone4est[1][:, 4, :]
    @test all(isapprox.(true_o_coords, est_o_coords, atol=0.2))

    """
    # The actual distribution of distances between the true and estimated oxygen atom positions:
 
    julia> histogram(mapslices(norm, true_o_coords .- est_o_coords, dims=1), nbins=9)
                 ┌                                        ┐ 
    [0.0 , 0.02) ┤████████████████████▌ 20                  
    [0.02, 0.04) ┤████████████████████████████████████  35  
    [0.04, 0.06) ┤██████████████████████████████████▉ 34    
    [0.06, 0.08) ┤█████████████████████▋ 21                 
    [0.08, 0.1 ) ┤█████████████████▌ 17                     
    [0.1 , 0.12) ┤██████████▍ 10                            
    [0.12, 0.14) ┤████████▎ 8                               
    [0.14, 0.16) ┤████▎ 4                                   
    [0.16, 0.18) ┤██▏ 2                                     
                 └                                        ┘

    # The estimate is on average ~0.05 angstroms away from the original PDB coords.
    """

end