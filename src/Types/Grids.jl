"""
    TransverseGrid(x, y)

A uniform transverse sampling grid for field representation.

# Arguments
- `x`: x-axis sample points [m], supplied as a range or vector
- `y`: y-axis sample points [m], supplied as a range or vector

# Fields
- `x`, `y`: coordinate vectors [m]
- `Nx`, `Ny`: number of sample points
- `dx`, `dy`: sample spacing [m]

# Examples
```jldoctest
julia> grid = TransverseGrid(
  range(-1e-3, 1e-3, length=256),
  range(-1e-3, 1e-3, length=256)
)
```
"""
struct TransverseGrid <: AbstractGrid
  x :: Vector{Float64}
  y :: Vector{Float64}
  Nx :: Int
  Ny :: Int
  dx :: Float64
  dy :: Float64

  function TransverseGrid(x, y)
    x = collect(x)
    y = collect(y)
    new(x, y, length(x), length(y), x[2]-x[1], y[2]-y[1])
  end
end
