module dihedrals
using LinearAlgebra
using Flux:batched_mul
using StatsBase
using Rotations

export dihedrals2xyz, dihedrals2xyz_exact, get_dihedrals, get_ks_ls_dihs, idealize_lengths_angles
bond_lengths, bond_angles = (1.4390292125539965, 2.0325377936819717) # Reasonable value for bond lengths and bond angles to be idealized to. (mean of approx 160k PDB values)

plotsh(x) = reshape(x, 3, :)

"""
Accepts a protein chain and returns the vectors and lengths connecting each pair of adjacent atoms. 
"""
function get_bond_vecs(protxyz)
    prot_flattened = plotsh(protxyz)
    len_prot = size(prot_flattened,2)
    pairwise_vecs = prot_flattened[:,2:len_prot] .- prot_flattened[:,1:len_prot-1]  
    lengths = sqrt.(pairwise_vecs[1,:].^2 .+ pairwise_vecs[2,:].^2 .+ pairwise_vecs[3,:].^2)
    return pairwise_vecs, lengths
end

"""
Turns a list of vectors into a set of points starting at the origin, where p1 = 0, p2 = v1 + p1, p3 = v2 + p2, etc.
"""
function coords_from_vecs(vecs)
    num_atoms = size(vecs,2) + 1
    coords = zeros(3,num_atoms)
    for i in 2:num_atoms
        coords[:,i] = coords[:,i-1] .+ vecs[:,i-1]
    end
    return coords
end

"""
Getting dihedral angles from vectors stored as a 3xN matrix and their lengths.
"""
function get_dihedrals(vecs, lengths)
    len_prot = size(vecs,2) +1

    crosses  = stack(cross.(eachcol(vecs[:,1:len_prot-2]), eachcol(vecs[:,2:len_prot-1])),dims=2)
    lengths_cross = reshape(sqrt.(crosses[1,:].^2 .+ crosses[2,:].^2 .+ crosses[3,:].^2),1,:)
    normalized_crosses = crosses ./ lengths_cross
    cos_theta = dot.(eachcol(normalized_crosses[:,1:len_prot-3]), eachcol(normalized_crosses[:,2:len_prot-2]))

    cross_crosses = stack(cross.(eachcol(normalized_crosses[:,1:len_prot-3]), eachcol(normalized_crosses[:,2:len_prot-2])),dims=2)
    normalized_vecs = (vecs ./ reshape(lengths,1,:))[:,2:len_prot-2]

    sin_theta = dot.(eachcol(cross_crosses), eachcol(normalized_vecs))
    
    thetas = atan.(sin_theta, cos_theta)
    return thetas
end

"""
Get dihedrals from protxyz only. 
"""
get_dihedrals(protxyz) = get_dihedrals(get_bond_vecs(protxyz)...)

"""
Law of cosines to find the third length of a triangle given two lengths and their intermediate angle. 
"""
function get_triangle(l1, l2, γ)
    l3 = sqrt(l1^2 + l2^2 - 2*l1*l2*cos(γ))
    return l3
end

"""
Fix the bond lengths of respective bonds to NCa, CaC, and CN, and return the new coordinates. 
"""
function fixed_bond_lengths(protxyz, NCa, CaC, CN)
    vecs, lengths = get_bond_vecs(protxyz)
    Nvecs = vecs ./ reshape(lengths,1,:) 
    new_coords = reshape(zeros(size(protxyz)),3,:)
    new_coords[:,1] = protxyz[:,1,1]
    for i in axes(Nvecs,2)
        if i % 3 == 1 
            new_coords[:,i+1] = new_coords[:,i] .+ NCa*Nvecs[:,i]
        elseif i % 3 == 2
            new_coords[:,i+1] = new_coords[:,i] .+ CaC*Nvecs[:,i]
        else
            new_coords[:,i+1] = new_coords[:,i] .+ CN*Nvecs[:,i]
        end
    end
    return reshape(new_coords,3,3,:)
end

"""
Idealizes the bond lengths and angles of protxyzs while maintaining the same overall structure. Protxyzs can be a single 3x3xN matrix or a vector of 3x3xN matrices.
"""
function idealize_lengths_angles(protxyzs::AbstractVector, bond_lengths = bond_lengths, bond_angles = bond_angles)

    prots_prepped = []

    # Code is general as to allow for different bond lengths and angles for different bonds, however it is currently not used. 
    NCa, CaC, CN = bond_lengths, bond_lengths, bond_lengths 
    NCaC, CaCN, CNCa = bond_angles, bond_angles, bond_angles 

    # get_triangle uses the law of cosines to find the length of the third side of a triangle given two bond lengths and one bond angle.
    lNCaC = get_triangle(NCa, CaC, NCaC)
    lCaCN = get_triangle(CaC, CN, CaCN)
    lCNCa = get_triangle(CN, NCa, CNCa)


    for i in axes(protxyzs, 1)
        org = fixed_bond_lengths(protxyzs[i], NCa, CaC, CN)
        orgp = protxyzs[i]
        points = plotsh(org[:,:,:])

        # ks and ls are lengths used for fixing bond angles. 
        ks = repeat([CaC, CN, NCa], size(orgp,3))
        ls = repeat([lCaCN, lCNCa, lNCaC], size(orgp,3))
        p = reshape(fix_sequence_of_points(points,ks,ls),3,3,:)
        push!(prots_prepped, p)
    end
    return prots_prepped
end

function idealize_lengths_angles(protxyz::AbstractArray, bond_lengths = bond_lengths, bond_angles = bond_angles)
    return idealize_lengths_angles([protxyz], bond_lengths, bond_angles)[1]
end

"""
Fixes P3 to the correct bond angle given its bond length and skip length, while keeping it in the same plane defined by the initial P1-P2-P3. 
"""
function fix_bond_angle(P1, P2, P3, k, l) 
    v1 = P1 - P2
    v2 = P3 - P2
    n = cross(v1,v2)
    p = P2
    d = norm(P1 - P2)

    # circle of intersection center
    h = (k^2-l^2+d^2) / (2*d)
    c = P2 + (h/d) * (P1-P2)

    # circle of intersection radius
    r = sqrt(k^2 - h^2)

    # plane's point closest to circle
    t = dot(n, (p-c)) / dot(n, n)
    proj = c + t * n

    # distance from the circle's center to the projection
    dist_to_proj = norm(proj - c)

    d = sqrt(r^2 - dist_to_proj^2)  # distance from projection to intersection points
    dir = cross(n, (P1 - P2))  # direction from proj to the intersection points
    dir /=norm(dir)  # normalize the direction

    # 2 intersection points
    P_1st = proj + d*dir
    P_2nd = proj - d*dir

    # choose one with right chirality
    i = argmin([norm(P_1st - P3), norm(P_2nd - P3 )])
    return [P_1st, P_2nd][i]
end

"""
Fix a sequence of points to their bond angles determined by ks and ls (bond lenths and skip lengths).
"""
function fix_sequence_of_points(points, ks, ls)
    # Initial checks
    new_points = deepcopy(points)
    for i in 3:size(points,2)
        new_points[:,i] = fix_bond_angle(new_points[:,i-2], new_points[:,i-1], new_points[:,i], ks[i-2], ls[i-2])
    end

    return new_points
end

"""
Edits the dihedrals of the given protxyz to a list or 3xN matrix of new dihedrals. 
"""
function dihedrals_to_vecs_respect_bond_angles(protxyz, new_dihedrals::AbstractVecOrMat)
    # Start at v_1. 
    # Compute dihedral angle v_1 - v2 - v_3
    # Update v_3 so that v_1 - v_2 - v_3 has the new dihedral angle. 
    # Continue for v_2 - v_3 - v_4 and so on...
    new_dihedrals = reshape(new_dihedrals,:)
    vectahs, lens = get_bond_vecs(protxyz)
    for i in 1:size(vectahs,2)-2 # Seems like there should exist a smart way to vectorize this

        # In every iteration we consider v_{i-1} - v_i - v_{i+1}
        curr_vecs = vectahs[:,i:i+2]
        lengths = reshape(mapslices(norm, curr_vecs, dims = 1),1,:)

        # Get the rotational update needed, in radians and about the v_i axis
        theta = mod.(new_dihedrals[i] .- get_dihedrals(curr_vecs,lengths) .+ π, 2 * π) .- π

        # Compute the rotation matrix corresponding to that rotation
        current_rotation = AngleAxis(theta..., curr_vecs[:,2]...)
        
        # Apply the rotation. 
        vectahs[:,i+2:end] = batched_mul(current_rotation , vectahs[:,i+2:end])
    end
    return vectahs, lens
end

""" 
Takes an array or a 3xN matrix of dihedrals and a starting residue and returns the xyz coordinates determined by the dihedrals, bond lengths and bond angles. 
"""
function dihedrals2xyz(dihs::AbstractVecOrMat, start_res::AbstractMatrix; bond_lengths = bond_lengths, bond_angles = bond_angles)
    dihs3xL = reshape(dihs,3,:) # note there are N-1 sets of three dihedrals for N residues. 
    reshape(start_res,3,3)
    init_points = cat(start_res, randn(3,3,size(dihs3xL,2)), dims = 3) 
    st = idealize_lengths_angles([init_points],bond_lengths,bond_angles)[1]
    
    return reshape(coords_from_vecs(dihedrals_to_vecs_respect_bond_angles(st, dihs3xL)[1]) .+ start_res[:,1],3,3,:)
end

"""
Returns the bond lengths between adjacent atoms (ks) and the skip lengths (ls).
"""
function get_ks_ls(protxyz)
    pxyz = reshape(protxyz, 3, :)
    _, ks = get_bond_vecs(pxyz)
    ls = []
    for i in 1:size(pxyz,2)-2
        push!(ls, norm(pxyz[:,i] .- pxyz[:,i+2]))
    end
    return ks, ls
end 

"""
Gets bond lengths, skip lengths, and dihedrals from protxyz. 
"""
function get_ks_ls_dihs(protxyz)
    ks, ls = get_ks_ls(protxyz)
    dihs = get_dihedrals(protxyz)
    return ks, ls, dihs
end

"""
Fix all bond lengths in protxyz to the bond lengths given in l. 
"""
function fix_bond_lengths_sequence(protxyz, ks)
    vecs, lengths = get_bond_vecs(protxyz)
    Nvecs = vecs ./ reshape(lengths,1,:) 
    new_coords = reshape(zeros(size(protxyz)),3,:)
    new_coords[:,1] = protxyz[:,1,1]
    for i in axes(Nvecs,2)
        ks[i]
        new_coords[:,i+1] = new_coords[:,i] .+ ks[i]*Nvecs[:,i]
    end
    return new_coords
end
""" 
Fixes the protxyzs to the exact sequence of lengths and angles (angles parametrized by skip lengths "ls"), where ks = [NCa1, CaC1, CN1,...] and ls = [NC, CaN, CCa,...].
"""
function fix_sequence_lengths_angles(protxyz, ks, ls)
    org = fix_bond_lengths_sequence(protxyz, ks)
    points = plotsh(org[:,:,:])
    p = reshape(fix_sequence_of_points(points,ks[1:end],ls),3,3,:)
    return p
end

"""
Maps dihedrals, adjacent bond lengths, skips lengths, and a starting residue to xyz coordinates. 
"""
function dihedrals2xyz_exact(dihs::AbstractVecOrMat, start_res::AbstractMatrix, ks::AbstractVector, ls::AbstractVector)
    dihs3xL = reshape(dihs,3,:) 
    init_points = cat(start_res, randn(3,3,size(dihs3xL,2)), dims = 3)
    st = fix_sequence_lengths_angles(init_points,ks,ls)
    return reshape(coords_from_vecs(dihedrals_to_vecs_respect_bond_angles(st, dihs3xL)[1]) .+ start_res[:,1],3,3,:)
end

end