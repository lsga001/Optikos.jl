using CairoMakie

@recipe FieldPlot (field::ScalarField,) begin
  "Defines the units of the plot"
  units = :mm
end

function Makie.plot!(fp::FieldPlot{<:Tuple{<:ScalarField}})
  input_nodes = [:field, :units]
  output_nodes = [:x, :y, :intensity]
  map!(fp.attributes, input_nodes, output_nodes) do field, units
    scale = units == :mm ? 1e3 :
            units == :μm ? 1e6 : 
            units == :nm ? 1e9 : 1.0
    x = @. field.grid.x * scale
    y = @. field.grid.y * scale
    I = intensity(field)
    return x, y, I
  end

  heatmap!(fp, fp.x, fp.y, fp.intensity)

  return fp 
end

Makie.convert_arguments(::Type{<:Makie.Plot{fieldplot}}, field::ScalarField) = (field,)
Makie.convert_arguments(::Type{<:Makie.Plot{fieldplot}}, grid::TransverseGrid{<:Real}, beam::AbstractBeam) = (evaluate(grid, beam),)

@recipe PhasePlot (field::ScalarField,) begin
    "Defines the units of the plot"
    units    = :mm
    "Colormap for phase display — should be cyclic"
    colormap = :hsv
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
      φ = angle.(field.U)
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
