"""
    TransverseGrid(x::Vector{Float64}[, y::Vector{Float64}])

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

julia> grid = TransverseGrid(x);
```
"""
struct TransverseGrid{T<:Real} <: AbstractGrid
  x :: Vector{T}
  y :: Vector{T}
  Nx :: Int
  Ny :: Int
  dx :: T
  dy :: T

  function TransverseGrid(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Real}
    x = collect(x)
    y = collect(y)
    new{T}(x, y, length(x), length(y), x[2]-x[1], y[2]-y[1])
  end
  
end

function TransverseGrid(x::AbstractVector{<:Real})
  return TransverseGrid(x, x)
end
