"""
    TransverseGrid(x::Vector{Float64}, y::Vector{Float64})

A uniform transverse sampling grid for field representation.

# Arguments
- `x`: x-axis sample points [m]
- `y`: y-axis sample points [m]

# Fields
- `x`, `y`: coordinate vectors [m] (Vector{Float64})
- `Nx`, `Ny`: number of sample points (Int)
- `dx`, `dy`: sample spacing [m] (Float64)

# Examples

```jldoctest
julia> x = range(-1e-3, 1e-3, length=256);

julia> y = range(-1e-3, 1e-3, length=256);

julia> grid = TransverseGrid(x, y);
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
