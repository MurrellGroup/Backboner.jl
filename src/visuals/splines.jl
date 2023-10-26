import Dierckx

function spline3d(points_matrix::AbstractMatrix{<:Real}; m::Integer=10, k::Integer=3)
    N = size(points_matrix, 2)
    spl = Dierckx.ParametricSpline(1:N, points_matrix, k=k)
    linrange = Dierckx.LinRange(1, N, m*N)
    points_matrix_fine = Dierckx.evaluate(spl, linrange)
    return points_matrix_fine
end

function spline3d(points::AbstractVector{<:AbstractVector{<:Real}}; kwargs...)
    points_matrix = hcat(points...)
    return spline3d(points_matrix; kwargs...)
end
