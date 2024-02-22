export is_knotted

"""
This file contains an implementation of the KnotFind algorithm, described by Khatib, Weirauch,
and Rohl in 2006, in the paper "Rapid knot detection and application to protein structure prediction"

Essentially:
 1. Take line segments formed by alphacarbons of a protein backbone.
 2. For an individual triple, i-1,i,i+1, if no line segments j,j+1 connecting consecutive Cα atoms
cross through the triangle, then Cα i is removed from the chain.
 3. If any line segment intersects the triangle, then no atom is removed
and the algorithm proceeds to the *next shortest*.
 4. After any Cα is removed, the algorithm returns to the triple with the shortest distance.
This is repeated until the last triple is reached and simplified, if possible.
 5. When it terminates, the protein contains no knots if there's only one segment left.
"""

function line_segment_intersects_triangle(
    p1::AbstractVector{T}, p2::AbstractVector{T}, a::AbstractVector{T}, b::AbstractVector{T}, c::AbstractVector{T},
) where T <: Real
    epsilon = eps(T)
    segment_origin = p1
    segment_vector = p2 - p1
    edge1 = b - a
    edge2 = c - a
    segment_cross_e2 = cross(segment_vector, edge2)
    det = dot(edge1, segment_cross_e2)
    abs(det) < epsilon && return false    # This segment is parallel to this triangle.
    inv_det = 1.0 / det
    s = segment_origin - a
    u = inv_det * dot(s, segment_cross_e2)
    (u < 0 || u > 1) && return false
    s_cross_e1 = cross(s, edge1)
    v = inv_det * dot(segment_vector, s_cross_e1)
    (v < 0 || u + v > 1) && return false
    t = inv_det * dot(edge2, s_cross_e1)
    return epsilon < t <= 1 # segment intersection
end

function splice_atom!(bond_vectors::Vector{Vector{T}}, bond_angles::Vector{T}, i::Integer) where T <: Real
    bond_vectors[i] = bond_vectors[i] + bond_vectors[i+1]
    splice!(bond_vectors, i+1)
    i > 1 && (bond_angles[i-1] = calculate_bond_angle(bond_vectors[i-1], bond_vectors[i]))
    i < length(bond_angles) && (bond_angles[i+1] = calculate_bond_angle(bond_vectors[i], bond_vectors[i+1]))
    splice!(bond_angles, i)
    return nothing
end

function get_bond_triangle(
    bond_vectors::Vector{Vector{T}}, i::Integer
) where T <: Real
    return zeros(T, 3), bond_vectors[i], bond_vectors[i] + bond_vectors[i+1]
end

function check_intersection(bond_vectors::Vector{Vector{T}}, i::Integer) where T <: Real
    a, b, c = get_bond_triangle(bond_vectors, i)
    
    p1 = c
    for j in i+2:length(bond_vectors)
        p2 = p1 + bond_vectors[j]
        line_segment_intersects_triangle(p1, p2, a, b, c) && return true
        p1 = p2
    end

    p1 = a
    for j in i-1:-1:1
        p2 = p1 - bond_vectors[j]
        line_segment_intersects_triangle(p1, p2, a, b, c) && return true
        p1 = p2
    end

    return false
end

function is_knotted(backbone::Backbone{T}) where T <: Real
    length(backbone) < 3 && return false

    bond_vectors = Vector.(eachcol(get_bond_vectors(backbone)))
    bond_angles = get_bond_angles(backbone)

    while !isempty(bond_angles)
        order = sortperm(bond_angles)
        has_spliced = false
        for i in order
            if !check_intersection(bond_vectors, i)
                splice_atom!(bond_vectors, bond_angles, i)
                has_spliced = true
                break
            end
        end
        !has_spliced && break
    end

    return !isempty(bond_angles)
end
