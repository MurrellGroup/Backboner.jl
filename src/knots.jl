"""
This file contains an implementation of the KnotFind algorithm, described by Khatib, Weirauch,
and Rohl in 2006, in the paper "Rapid knot detection and application to protein structure prediction"

 1. Take all alphacarbons of a protein backbone.
 2. Sort them using some metric, e.g. distance or area of the triangle formed by three consecutive Cα atoms.
 3. For an individual triple, i-1,i,i+1, if no line segments j,j+1 connecting consecutive Cα atoms
cross through the triangle, then Cα i is removed from the chain.
 4. If any line segment intersects the triangle, then no atom is removed
and the algorithm proceeds to the triangle with the next smallest metric value.
 4. After any Cα is removed, the algorithm returns to the triple with the shortest distance, just in case some segment no longer intersects,
This is repeated until the last triple is reached and simplified, if possible.
 5. When it terminates, the protein contains no knots if there's only one segment left.
 6. Run again with a different metric if the chain is not simplified to make sure it wasn't a false positive.
"""

using StaticArrays: SVector # leads to an insane speed-up

triangle_distance(a, _, c) = norm(c - a)
triangle_area(a, b, c) = norm(cross(b - a, c - a)) / 2

triangle_values(metric::Function, points::Vector{<:AbstractVector{T}}) where T<:Real = T[metric(points[i:i+2]...) for i in 1:length(points)-2]

function line_segment_intersects_triangle(segment_start, segment_end, a, b, c)
    T = eltype(segment_start)
    segment_vector = segment_end - segment_start
    edge1, edge2 = b - a, c - a
    segment_cross_e2 = cross(segment_vector, edge2)
    det = dot(edge1, segment_cross_e2)
    abs(det) < eps(T) && return false # segment is parallel to triangle
    inv_det = 1.0 / det
    s = segment_start - a
    u = inv_det * dot(s, segment_cross_e2)
    (u < 0 || u > 1) && return false
    s_cross_e1 = cross(s, edge1)
    v = inv_det * dot(segment_vector, s_cross_e1)
    (v < 0 || u + v > 1) && return false
    t = inv_det * dot(edge2, s_cross_e1)
    return eps(T) < t <= 1 # segment intersection
end

function check_intersection(points::Vector{<:AbstractVector{T}}, i::Int) where T<:Real
    a, b, c = points[i-1], points[i], points[i+1] 
    for j in 1:length(points)-1
        i-2 <= j <= i+1 && continue
        p1, p2 = points[j], points[j+1]
        line_segment_intersects_triangle(p1, p2, a, b, c) && return true
    end
    return false
end

# for removing an atom from a backbone, and updating adjacent triangles
function remove_atom!(points::Vector{<:AbstractVector{T}}, i::Int, metric_values::Vector{T}, metric::Function) where T<:Real
    m = length(metric_values)
    triangle_index = i - 1
    triangle_index > 1 && (metric_values[triangle_index-1] = metric(points[i-2], points[i-1], points[i+1]))
    triangle_index < m && (metric_values[triangle_index+1] = metric(points[i-1], points[i+1], points[i+2]))
    deleteat!(points, i)
    deleteat!(metric_values, triangle_index)
    return nothing
end

function simplify!(points::Vector{<:AbstractVector{T}}, metric::Function) where T<:Real
    metric_values = triangle_values(metric, points) # Vector{T} of length n-2
    while !isempty(metric_values)
        order = sortperm(metric_values) # TODO: calculate once outside while loop, update in `remove_atom!`
        has_removed = false
        for triangle_index in order
            i = triangle_index + 1
            if !check_intersection(points, i)
                remove_atom!(points, i, metric_values, metric)
                has_removed = true
                break
            end
        end
        !has_removed && break # terminate if the chain can't be simplified further
    end
    return points
end

simplify(points::Vector{<:AbstractVector{T}}, metric::Function) where T<:Real = simplify!(deepcopy(points), metric)

"""
    is_knotted(backbone::Backbone)

Check if a backbone is knotted.
"""
function is_knotted(backbone::Backbone, metrics::Vector{Function} = [triangle_distance, triangle_area])
    length(backbone) < 2 && return false
    points = SVector{3}.(backbone) # convert to StaticArrays for 40x performance lmao
    for metric in metrics
        length(simplify(points, metric)) == 2 && return false
    end
    return true
end
