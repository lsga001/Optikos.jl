"""
    realize(source::DiskSchellModel, grid::TransverseGrid) -> ScalarField

Draw a single random realization from a `DiskSchellModel` source, evaluated at a 
transverse plane a distance `z` [m] from the source plane.
Returns a paraboloidal wave originating from a point drawn uniformly
at random from within a disk of radius `R`.

Each call is statistically independent.

Each sample is a paraboloidal wave:

```math
E(x,y) = \\frac{1}{z}\\exp(ikz)\\exp\\!\\left(ik\\frac{(x-x_0)^2+(y-y_0)^2}{2z}\\right)
```

where `(x₀, y₀)` is drawn uniformly from a disk of radius `R`.

# Arguments
- `source`: `DiskSchellModel` source descriptor
- `grid`  : transverse observation grid
- `z`     : propagation distance from source to observation plane [m]
"""
function sample(
  source::DiskSchellModel, 
  grid::TransverseGrid{<:Real}, 
  z::Float64
) :: Tuple{ScalarField, Tuple{Float64, Float64}}

  x0, y0 = _sample_disk(source.R)
  point = (x0, y0)
 
  λ = source.λ
  beam = ParaboloidalBeam(λ, z, center=(x0, y0))
  field = evaluate(grid, beam)

  #return ScalarField(E, grid, source.λ)
  return field, point
end

# Internal helper — uniform sampling within a disk via rejection sampling
function _sample_disk(R::Float64) :: Tuple{Float64, Float64}
  r = R * sqrt(rand())
  φ = 2π * rand()
  return r * cos(φ), r * sin(φ)
end
