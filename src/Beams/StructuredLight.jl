"""
    evaluate(beam::GaussianBeam, grid::TransverseGrid) -> ScalarField

Sample a Gaussian beam onto a transverse grid at the plane `z = z0`.

The complex field amplitude at the waist (`z = 0`) is:

```math
E(x, y) = \\exp\\!\\left(-\\frac{x^2 + y^2}{w_0^2}\\right)
```

For `z ≠ 0`, the full Gaussian beam expression including the beam radius
`w(z)`, wavefront curvature `R(z)`, and Gouy phase `ψ(z)` is used:

```math
E(x, y, z) = \\frac{w_0}{w(z)}
\\exp\\!\\left(-\\frac{r^2}{w(z)^2}\\right)
\\exp\\!\\left(-i\\left(kz + \\frac{kr^2}{2R(z)} - \\psi(z)\\right)\\right)
```

where:
- ``w(z) = w_0\\sqrt{1 + (z/z_R)^2}``
- ``R(z) = z\\left(1 + (z_R/z)^2\\right)``
- ``\\psi(z) = \\arctan(z/z_R)``
- ``z_R = \\pi w_0^2 n / \\lambda`` is the Rayleigh range

# Arguments
- `beam`: [`GaussianBeam`](@ref) descriptor
- `grid`: [`TransverseGrid`](@ref) on which to sample

# Returns
A [`ScalarField`](@ref) with unit peak amplitude at the waist.

# Examples
```jldoctest
julia> grid  = TransverseGrid(range(-1e-3, 1e-3, 64), range(-1e-3, 1e-3, 64));

julia> beam  = GaussianBeam(200e-6, 632.8e-9);

julia> field = evaluate(beam, grid);

julia> size(field.E)
(64, 64)

julia> abs(field.E[32, 32]) < 1.0   # off-center, less than peak
true
```

# References
Saleh & Teich, *Fundamentals of Photonics*, 3rd ed., §3.1
"""
function evaluate(beam::GaussianBeam, grid::TransverseGrid) :: ScalarField
  (; w0, λ, z0, n_index) = beam

  k = 2π * n_index / λ
  zR = π * w0^2 * n_index / λ # Rayleigh range
  z = z0

  X = grid.x' .* ones(grid.Ny)
  Y = ones(grid.Nx)' .* grid.y
  r² = X.^2 .+ Y.^2

  U = if iszero(z)
    @. exp(-r² / w0^2) + 0im
  else
    wz = w0 * sqrt(1 + (z/zR)^2 ) # beam waist at z
    Rz = z * (1 + (zR / z)^2) # wavefront radius at z
    ψz = atan(z / zR) # Gouy phase
    @. (w0 / wz) *
      exp(-r² / wz^2) *
      exp(-1im * (k*z + k*r² / (2*Rz) - ψz))
  end

  return ScalarField(U, grid, λ)
end
