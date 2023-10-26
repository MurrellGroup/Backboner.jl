include("splines.jl")


function plot_spline!(backbone::BackboneAndOxygen; linewidth=5, color=:blue)
    CA_coords = alphacarbon_coord_matrix(backbone)
    CA_coords_fine = spline3d(CA_coords)
    plot!(eachrow(CA_coords_fine)..., linewidth=linewidth, color=color)
end

function plot_sheet!(backbone::BackboneAndOxygen; linewidth=15, color=:brown)
    CA_coords = alphacarbon_coord_matrix(backbone)
    CA_coords_fine = spline3d(CA_coords)
    plot!(eachrow(CA_coords_fine)..., linewidth=linewidth, color=color)
end


function plot_segment!(segment::Segment{Loop}, backbone::BackboneAndOxygen)
    plot_spline!(backbone[segment.range], linewidth=6, color=:tan)
end

function plot_segment!(segment::Segment{Helix}, backbone::BackboneAndOxygen)
    plot_spline!(backbone[segment.range], linewidth=15, color=:brown)
end

function plot_segment!(segment::Segment{Strand}, backbone::BackboneAndOxygen)
    plot_sheet!(backbone[segment.range], linewidth=15, color=:green)
end


function plot_segments!(segments::AbstractVector{<:Segment}, backbone::BackboneAndOxygen)
    for segment in segments
        plot_segment!(segment, backbone)
    end
end