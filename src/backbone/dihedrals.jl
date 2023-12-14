using LinearAlgebra
using NNlib: batched_mul
using Rotations

export dihedrals2xyz, dihedrals2xyz_exact, get_dihedrals, get_ks_ls, get_ks_ls_dihs, idealize_lengths_angles, new_frame_dihedrals

# Mean value of bond lengths and bond angles, from approx. 160k PDB values
const MEAN_BOND_LENGTH = 1.439029
const MEAN_BOND_ANGLE = 2.032538

#=
TODO: think about how to better integrate this with the Backbone type.
could also probably distill the set of functions to be more concise and general,
e.g. 
=#

"""
Returns the vectors and lengths connecting each pair of adjacent atoms in the backbone
"""
function bonds_vecs_and_lens(backbone::Backbone{N}) where N
    @assert N >= 3 "backbone needs at least the N, Cα and C atoms to calculate bond vectors"
    bond_vectors = backbone_bond_vectors(backbone)
    lengths = reshape(mapslices(norm, bond_vectors, dims=1), :)
    return bond_vectors, lengths
end


"""
Turns a list of vectors into a set of points starting at the origin, where p1 = 0, p2 = v1 + p1, p3 = v2 + p2, etc.
"""
function coords_from_vecs(vecs)
    num_atoms = size(vecs, 2) + 1
    coords = zeros(3, num_atoms)
    for i in 2:num_atoms
        coords[:, i] = coords[:, i-1] .+ vecs[:, i-1]
    end
    return coords
end

"""
    get_dihedrals(vecs::AbstractMatrix, lengths::AbstractVector)

Getting dihedral angles from vectors stored as a 3xN matrix and their lengths.
"""
function get_dihedrals(vecs::AbstractMatrix, lengths::AbstractVector)
    len_prot = size(vecs, 2) + 1

    crosses = stack(cross.(eachcol(vecs[:, 1:len_prot-2]), eachcol(vecs[:, 2:len_prot-1])), dims=2)
    lengths_cross = reshape(sqrt.(crosses[1, :].^2 .+ crosses[2, :].^2 .+ crosses[3, :].^2), 1, :)
    normalized_crosses = crosses ./ lengths_cross
    cos_theta = dot.(eachcol(normalized_crosses[:, 1:len_prot-3]), eachcol(normalized_crosses[:, 2:len_prot-2]))

    cross_crosses = stack(cross.(eachcol(normalized_crosses[:, 1:len_prot-3]), eachcol(normalized_crosses[:, 2:len_prot-2])), dims=2)
    normalized_vecs = (vecs ./ reshape(lengths, 1, :))[:, 2:len_prot-2]

    sin_theta = dot.(eachcol(cross_crosses), eachcol(normalized_vecs))
    
    thetas = atan.(sin_theta, cos_theta)
    return thetas
end

"""
    get_dihedrals(backbone::Backbone)

Getting dihedral angles from a Backbone.
"""
@inline get_dihedrals(backbone::Backbone) = get_dihedrals(bonds_vecs_and_lens(backbone)...)
@inline get_dihedrals(coords::AbstractArray{T, 3}) where T = get_dihedrals(Backbone(coords))

"""
Fix the bond lengths of respective bonds to NCa, CaC, and CN, and return the new coordinates.
"""
function fixed_bond_lengths(coords::AbstractArray{T, 3}, NCa, CaC, CN) where T
    vecs, lengths = bonds_vecs_and_lens(Backbone(coords))
    Nvecs = vecs ./ reshape(lengths, 1, :) 
    new_coords = reshape(zeros(size(coords)), 3, :)
    new_coords[:, 1] = coords[:, 1, 1]
    for i in axes(Nvecs, 2)
        if i % 3 == 1 
            new_coords[:, i+1] = new_coords[:, i] .+ NCa*Nvecs[:, i]
        elseif i % 3 == 2
            new_coords[:, i+1] = new_coords[:, i] .+ CaC*Nvecs[:, i]
        else
            new_coords[:, i+1] = new_coords[:, i] .+ CN*Nvecs[:, i]
        end
    end
    return reshape(new_coords, 3, 3, :)
end

# Finds the length of the third side of a triangle given two lengths and their intermediate angle.
@inline law_of_cosines(a, b, C) = sqrt(a^2 + b^2 - 2*a*b*cos(C))

"""
    idealize_lengths_angles(coords_vector::AbstractVector{<:AbstractArray{T, 3}}; bond_lengths=MEAN_BOND_LENGTH, bond_angles=MEAN_BOND_ANGLE) where T

Idealizes the bond lengths and angles of coords_vector while maintaining the same overall structure. coords_vector can be a single 3x3xN matrix or a vector of 3x3xN matrices.
"""
function idealize_lengths_angles(
    coords_vector::AbstractVector{<:AbstractArray{T, 3}};
    bond_lengths=MEAN_BOND_LENGTH,
    bond_angles=MEAN_BOND_ANGLE,
) where T
    # Code is general as to allow for different bond lengths and angles for different bonds, however it is currently not used. 
    NCa = CaC = CN = bond_lengths 
    NCaC = CaCN = CNCa = bond_angles 

    lNCaC = law_of_cosines(NCa, CaC, NCaC)
    lCaCN = law_of_cosines(CaC, CN, CaCN)
    lCNCa = law_of_cosines(CN, NCa, CNCa)

    prots_prepped = Array{T, 3}[]

    for i in axes(coords_vector, 1)
        org = fixed_bond_lengths(coords_vector[i], NCa, CaC, CN)
        orgp = coords_vector[i]
        points = reshape(org[:, :, :], :)

        # ks and ls are lengths used for fixing bond angles. 
        ks = repeat([CaC, CN, NCa], size(orgp, 3))
        ls = repeat([lCaCN, lCNCa, lNCaC], size(orgp, 3))
        p = reshape(fix_sequence_of_points(points, ks, ls), 3, 3, :)
        push!(prots_prepped, p)
    end
    return prots_prepped
end

function idealize_lengths_angles(
    coords::AbstractArray{T, 3};
    bond_lengths=MEAN_BOND_LENGTH,
    bond_angles=MEAN_BOND_ANGLE
) where T
    return idealize_lengths_angles([coords]; bond_lengths=bond_lengths, bond_angles=bond_angles)[1]
end

"""
Fixes P3 to the correct bond angle given its bond length and skip length, while keeping it in the same plane defined by the initial P1-P2-P3.
"""
function fix_bond_angle(P1, P2, P3, k, l) 
    v1 = P1 - P2
    v2 = P3 - P2
    n = cross(v1, v2)
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
    dir /= norm(dir)  # normalize the direction

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
    for i in 3:size(points, 2)
        new_points[:, i] = fix_bond_angle(new_points[:, i-2], new_points[:, i-1], new_points[:, i], ks[i-2], ls[i-2])
    end

    return new_points
end

"""
Edits the dihedrals of the given coords to a list or 3xN matrix of new dihedrals.
"""
function dihedrals_to_vecs_respect_bond_angles(
    coords::AbstractArray{T, 3},
    new_dihedrals::AbstractVecOrMat
) where T
    # Start at v_1. 
    # Compute dihedral angle v_1 - v2 - v_3
    # Update v_3 so that v_1 - v_2 - v_3 has the new dihedral angle. 
    # Continue for v_2 - v_3 - v_4 and so on...
    new_dihedrals = reshape(new_dihedrals, :)
    vectahs, lens = bonds_vecs_and_lens(Backbone(coords))
    for i in 1:size(vectahs, 2)-2 # Seems like there should exist a smart way to vectorize this

        # In every iteration we consider v_{i-1} - v_i - v_{i+1}
        curr_vecs = vectahs[:, i:i+2]
        lengths = reshape(mapslices(norm, curr_vecs, dims=1), :)

        # Get the rotational update needed, in radians and about the v_i axis
        theta = mod.(new_dihedrals[i] .- get_dihedrals(curr_vecs, lengths) .+ π, 2π) .- π

        # Compute the rotation matrix corresponding to that rotation
        current_rotation = AngleAxis(theta..., curr_vecs[:, 2]...)

        # Apply the rotation. 
        vectahs[:, i+2:end] = batched_mul(current_rotation, vectahs[:, i+2:end])
    end
    return vectahs, lens
end

""" 
    dihedrals2xyz(dihedrals::AbstractVecOrMat, start_res::AbstractMatrix; bond_lengths=MEAN_BOND_LENGTH, bond_angles=MEAN_BOND_ANGLE)

Takes an array or a 3xN matrix of dihedrals and a starting residue and returns the xyz coordinates determined by the dihedrals, bond lengths and bond angles. 
"""
function dihedrals2xyz(
    dihedrals::AbstractVecOrMat,
    start_res::AbstractMatrix;
    bond_lengths=MEAN_BOND_LENGTH,
    bond_angles=MEAN_BOND_ANGLE,
)
    dihedrals3xL = reshape(dihedrals, 3, :) # note there are N-1 sets of three dihedrals for N residues.
    init_points = cat(start_res, randn(3, 3, size(dihedrals3xL, 2)), dims = 3) 
    st = idealize_lengths_angles([init_points], bond_lengths=bond_lengths, bond_angles=bond_angles)[1]
    return reshape(coords_from_vecs(dihedrals_to_vecs_respect_bond_angles(st, dihedrals3xL)[1]) .+ start_res[:, 1], 3, 3, :)
end

"""
    get_ks_ls(coords::AbstractArray{T, 3}) where T

Returns the bond lengths between adjacent atoms (ks) and the skip lengths (ls).
"""
function get_ks_ls(coords::AbstractArray{T, 3}) where T
    _, ks = bonds_vecs_and_lens(Backbone(coords))
    flat_coords = reshape(coords, 3, :)
    ls = [norm(@view(flat_coords[:, i]) .- @view(flat_coords[:, i+2])) for i in 1:size(flat_coords, 2)-2]
    return ks, ls
end 

"""
    get_ks_ls_dihs(coords::AbstractArray{T, 3}) where T

Gets bond lengths, skip lengths, and dihedrals from coords. 
"""
function get_ks_ls_dihs(coords::AbstractArray{T, 3}) where T
    ks, ls = get_ks_ls(coords)
    dihedrals = get_dihedrals(coords)
    return ks, ls, dihedrals
end

"""
Fix all bond lengths in coords to the bond lengths given in l.
"""
function fix_bond_lengths_sequence(coords::AbstractArray{T, 3}, ks::AbstractVector) where T
    vecs, lengths = bonds_vecs_and_lens(Backbone(coords))
    Nvecs = vecs ./ reshape(lengths, 1, :) 
    new_coords = reshape(zeros(size(coords)), 3, :)
    new_coords[:, 1] = coords[:, 1, 1]
    for i in axes(Nvecs, 2)
        new_coords[:, i+1] = new_coords[:, i] .+ ks[i]*Nvecs[:, i]
    end
    return new_coords
end
 
"""
Fixes the coords_vector to the exact sequence of lengths and angles (angles parametrized by skip lengths "ls"), where ks = [NCa1, CaC1, CN1, ...] and ls = [NC, CaN, CCa, ...].
"""
function fix_sequence_lengths_angles(coords::AbstractArray{T, 3}, ks::AbstractVector, ls::AbstractVector) where T
    org = fix_bond_lengths_sequence(coords, ks)
    points = reshape(org[:, :, :], :)
    p = reshape(fix_sequence_of_points(points, ks[1:end], ls), 3, 3, :)
    return p
end

"""
    dihedrals2xyz_exact(dihedrals::AbstractVecOrMat, start_res::AbstractMatrix, ks::AbstractVector, ls::AbstractVector)

Maps dihedrals, adjacent bond lengths, skips lengths, and a starting residue to xyz coordinates. 
"""
function dihedrals2xyz_exact(dihedrals::AbstractVecOrMat, start_res::AbstractMatrix, ks::AbstractVector, ls::AbstractVector)
    dihedrals3xL = reshape(dihedrals, 3, :) 
    init_points = cat(start_res, randn(3, 3, size(dihedrals3xL, 2)), dims = 3)
    st = fix_sequence_lengths_angles(init_points, ks, ls)
    return reshape(coords_from_vecs(dihedrals_to_vecs_respect_bond_angles(st, dihedrals3xL)[1]) .+ start_res[:, 1], 3, 3, :)
end

"""
    new_frame_dihedrals(frames_prev::AbstractArray{T, 4}, dihedrals::AbstractMatrix) where T

Takes protxyz and a singular set of dihedrals to place the next frame with idealized bond lengths and angles
"""
function new_frame_dihedrals(frames_prev::AbstractArray{T, 3}, dihedrals::AbstractMatrix) where T
    new_frame = dihedrals2xyz(dihedrals[:, end], frames_prev[:, :, end])
    pxyz = cat(frames_prev[:, :, :], new_frame[:, :, end], dims=3)
    return pxyz
end
