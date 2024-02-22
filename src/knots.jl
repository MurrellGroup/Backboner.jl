export is_knotted

#=
This file contains an implementation of the KnotFind algorithm, described by Khatib, Weirauch,
and Rohl in 2006, in the paper "Rapid knot detection and application to protein structure prediction"

 1. Take line segments formed by alphacarbons of a protein backbone.
 2. For an individual triple, i-1,i,i+1, if no line segments j,j+1 connecting consecutive Cα atoms
cross through the triangle, then Cα i is removed from the chain.
 3. If any line segment intersects the triangle, then no atom is removed
and the algorithm proceeds to the *next shortest*.
 4. After any Cα is removed, the algorithm returns to the triple with the shortest distance.
This is repeated until the last triple is reached and simplified, if possible.
 5. When it terminates, the protein contains no knots if there's only one segment left.
=#

@inline triangle_area(u::V, v::V) where {T <: Real, V <: AbstractVector{T}} = T(0.5) * norm(cross(u, v))
@inline triangle_area(a::V, b::V, c::V) where V <: AbstractVector{<:Real} = triangle_area(b - a, c - a)
triangle_areas(points::Vector{Vector{T}}) where T <: Real = [triangle_area(points[i:i+2]...) for i in 1:length(points)-2]

function line_segment_intersects_triangle(
    segment_start::AbstractVector{T}, segment_end::AbstractVector{T},
    a::AbstractVector{T}, b::AbstractVector{T}, c::AbstractVector{T},
) where T <: Real
    epsilon = eps(T)
    segment_vector = segment_end - segment_start
    edge1, edge2 = b - a, c - a
    segment_cross_e2 = cross(segment_vector, edge2)
    det = dot(edge1, segment_cross_e2)
    abs(det) < epsilon && return false # segment is parallel to this triangle.
    inv_det = 1.0 / det
    s = segment_start - a
    u = inv_det * dot(s, segment_cross_e2)
    (u < 0 || u > 1) && return false
    s_cross_e1 = cross(s, edge1)
    v = inv_det * dot(segment_vector, s_cross_e1)
    (v < 0 || u + v > 1) && return false
    t = inv_det * dot(edge2, s_cross_e1)
    return epsilon < t <= 1 # segment intersection
end

# for removing an atom from a backbone, making sure to update the surrounding points
function remove_atom!(points::Vector{Vector{T}}, areas::Vector{T}, i::Int) where T <: Real
    triangle_index = i - 1
    triangle_index > 1 && (areas[triangle_index-1] = triangle_area(points[i-2], points[i-1], points[i+1]))
    triangle_index < length(areas) && (areas[triangle_index+1] = triangle_area(points[i-1], points[i+1], points[i+2]))
    deleteat!(points, i)
    deleteat!(areas, triangle_index)
    return nothing
end

function check_intersection(points::Vector{Vector{T}}, i::Int) where T <: Real
    a, b, c = points[i-1], points[i], points[i+1] 
    for j in 1:length(points)-1
        i-1 <= j < i+1 && continue
        p1, p2 = points[j], points[j+1]
        line_segment_intersects_triangle(p1, p2, a, b, c) && return true
    end
    return false
end

function simplify!(points::Vector{Vector{T}}) where T <: Real
    areas = triangle_areas(points)
    while !isempty(areas)
        order = sortperm(areas) # TODO: calculate once, update in `remove_atom!`
        has_removed = false
        for triangle_index in order
            i = triangle_index + 1
            if !check_intersection(points, i)
                remove_atom!(points, areas, i)
                has_removed = true
                break
            end
        end
        !has_removed && break # halt if the chain couldn't be simplified further
    end
    return nothing
end

"""
    is_knotted(backbone::Backbone)

Check if a backbone is knotted.
"""
function is_knotted(backbone::Backbone{T}) where T <: Real
    points = Vector.(collect(eachcol(backbone.coords)))
    simplify!(points)
    return length(points) > 2
end
