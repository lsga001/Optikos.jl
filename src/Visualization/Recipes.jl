using CairoMakie

@recipe FieldPlot (field::ScalarField,) begin
  "Defines the units of the plot"
  units = :mm
  "Type of plot"
  plottype = :intensity
  colormap = :inferno
end

function Makie.plot!(fp::FieldPlot{<:Tuple{<:ScalarField}})
  input_nodes = [:field, :plottype, :units]
  output_nodes = [:x, :y, :intensity, :colorrange]
  map!(fp.attributes, input_nodes, output_nodes) do field, plottype, units
    scale = units == :mm ? 1e3 :
            units == :μm ? 1e6 : 
            units == :nm ? 1e9 : 1.0
    x = @. field.grid.x * scale
    y = @. field.grid.y * scale
    U = plottype == :real ? real(field)' :
        plottype == :imaginary ? imaginary(field)' : 
        plottype == :abs ? abs(field)' : intensity(field)'

    # Make it so that it doesn't fail for I_min≈I_max due to color
    U_min, U_max = extrema(U)
    colorrange = U_min ≈ U_max ? (0, 2*U_max) : (U_min, U_max)
    return x, y, U, colorrange
  end

  heatmap!(fp, fp.x, fp.y, fp.intensity, colorrange = fp.colorrange, colormap = fp.colormap)

  return fp 
end

Makie.convert_arguments(::Type{<:Makie.Plot{fieldplot}}, field::ScalarField) = (field,)
Makie.convert_arguments(::Type{<:Makie.Plot{fieldplot}}, grid::TransverseGrid{<:Real}, beam::AbstractBeam) = (evaluate(grid, beam),)

@recipe PhasePlot (field::ScalarField,) begin
    "Defines the units of the plot"
    units    = :mm
    "Colormap for phase display — should be cyclic"
    colormap = :twilight
end

function Makie.plot!(pp::PhasePlot{<:Tuple{<:ScalarField}})
  input_nodes = [:field, :units]
  output_nodes = [:x, :y, :phase]
  map!(pp.attributes, input_nodes, output_nodes) do field, units
      scale = units == :mm ? 1e3 :
              units == :μm ? 1e6 :
              units == :nm ? 1e9 : 1.0
      x = @. field.grid.x * scale
      y = @. field.grid.y * scale
      φ = angle.(field.U)'
      return x, y, φ
  end

  heatmap!(pp, pp.x, pp.y, pp.phase,
      colormap   = pp.colormap,
      colorrange = (-π, π),
  )
  return pp
end

Makie.convert_arguments(::Type{<:Makie.Plot{phaseplot}}, field::ScalarField) = (field,)
Makie.convert_arguments(::Type{<:Makie.Plot{phaseplot}}, grid::TransverseGrid{<:Real}, beam::AbstractBeam) = (evaluate(grid, beam),)
