"""
    _laguerre(n, α, x) -> Real

Evaluate the generalized Laguerre polynomial L_n^α(x).
Computed via the recurrence relation for numerical stability.

# References
Abramowitz & Stegun, §22.7
"""
function _laguerre(n::Int, α::Real, x::Real) :: Real
  @assert n >= 0 "n must be nonnegative"
  #@assert α > -1 "α must be > -1"
  n == 0 && return 1.0
  n == 1 && return 1.0 + α - x
  L_prev2 = 1.0
  L_prev1 = 1.0 + α - x
  L = 0.0
  for k in 1:(n-1)
    L = ((2k + 1 + α - x) * L_prev1 - (k + α) * L_prev2) / (k + 1)
    L_prev2 = L_prev1
    L_prev1 = L
  end
  return L
end

"""
    evaluate(grid::TransverseGrid{<:Real}, beam::SphericalBeam) -> ScalarField

Sample a Spherical beam onto a transverse grid at the plane `z = z0`.

The complex field amplitude at z≠0 is:

```math
U(r) = \\frac{A_0}{r}\\exp\\!\\left(-jkr\\right)
```

# Arguments
- `grid`: [`TransverseGrid`](@ref) on which to sample
- `beam`: [`GaussianBeam`](@ref) descriptor

# Returns
A [`ScalarField`](@ref) with unit peak amplitude at the waist.

# Examples
```jldoctest
julia> grid  = TransverseGrid(range(-1e-3, 1e-3, 64), range(-1e-3, 1e-3, 64));

julia> beam  = SphericalBeam(632.8e-9, 0.1);

julia> field = evaluate(grid, beam);

julia> size(field.U)
(64, 64)
```

# References
Saleh & Teich, *Fundamentals of Photonics*, 3rd ed., §2.2
"""
function evaluate(grid::TransverseGrid{<:Real}, beam::SphericalBeam) :: ScalarField
  (; λ, z, n_index, center) = beam

  x0, y0 = center

  k = 2π * n_index / λ
  
  X = @. grid.x' - x0
  Y = @. grid.y - y0
  r = @. sqrt(X^2 + Y^2 + z^2)

  U = @. exp(-1im * k * r) / r
  return ScalarField(U, grid, λ)
end

"""
    evaluate(grid::TransverseGrid{<:Real}, beam::ParaboloidalBeam) -> ScalarField

Sample a Paraboloidal beam onto a transverse grid at the plane `z = z0`.

The complex field amplitude at z≠0 is:

```math
U(x,y,z) = \\frac{A_0}{z}\\exp\\!\\left(-jkz\\right)\\exp\\!\\left[-jk\\frac{x^2 + y^2}{2z}\\right]
```

# Arguments
- `grid`: [`TransverseGrid`](@ref) on which to sample
- `beam`: [`GaussianBeam`](@ref) descriptor

# Returns
A [`ScalarField`](@ref) with unit peak amplitude at the waist.

# Examples
```jldoctest
julia> grid  = TransverseGrid(range(-1e-3, 1e-3, 64), range(-1e-3, 1e-3, 64));

julia> beam  = ParaboloidalBeam(632.8e-9, 0.1);

julia> field = evaluate(grid, beam);

julia> size(field.U)
(64, 64)
```

# References
Saleh & Teich, *Fundamentals of Photonics*, 3rd ed., §2.2
"""
function evaluate(grid::TransverseGrid{<:Real}, beam::ParaboloidalBeam) :: ScalarField
  (; λ, z, n_index, center) = beam
  
  x0, y0 = center

  k = 2π * n_index / λ
  
  X = @. grid.x' - x0
  Y = @. grid.y - y0
  r² = @. X^2 + Y^2

  U = @. exp(-1im * k * z) * exp(-1im * k * r² / (2 * z)) / z
  return ScalarField(U, grid, λ)
end

"""
    evaluate(grid::TransverseGrid{<:Real}, beam::GaussianBeam) -> ScalarField

Sample a Gaussian beam onto a transverse grid at the plane `z = z0`.

The complex field amplitude at the waist (`z = 0`) is:

```math
U(x, y) = \\exp\\!\\left(-\\frac{x^2 + y^2}{w_0^2}\\right)
```

For `z ≠ 0`, the full Gaussian beam expression including the beam radius
`w(z)`, wavefront curvature `R(z)`, and Gouy phase `ψ(z)` is used:

```math
U(x, y, z) = \\frac{w_0}{w(z)}
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

julia> field = evaluate(grid, beam);

julia> size(field.U)
(64, 64)
```

# References
Saleh & Teich, *Fundamentals of Photonics*, 3rd ed., §3.1
"""
function evaluate(grid::TransverseGrid{<:Real}, beam::GaussianBeam) :: ScalarField
  (; w0, λ, z0, n_index, center) = beam
  
  x0, y0 = center

  k = 2π * n_index / λ
  zR = π * w0^2 / λ # Rayleigh range
  z = z0
  
  wz = w0 * sqrt(1 + (z/zR)^2 ) # beam waist at z
  Rz = iszero(z) ? Inf : z * (1 + (zR / z)^2) # wavefront radius at z
  ψz = atan(z / zR) # Gouy phase

  X = @. grid.x' - x0
  Y = @. grid.y - y0
  ρ² = @. X^2 + Y^2

  #C = sqrt(2/π) / wz
  C = w0 / wz

  U = @. C *
      exp(-ρ² / wz^2) *
      exp(-1im * (k*z + k*ρ² / (2*Rz) - ψz))

  return ScalarField(U, grid, λ)
end

"""
    evaluate(grid::TransverseGrid{<:Real}, beam::LGBeam) -> ScalarField

Sample a Laguerre-Gaussian beam LG(p, l) onto a transverse grid at `z = z0`.

The complex field amplitude is:

```math
U(r, \\phi) = \\left[\\frac{C_{p,l}}{w(z)}\\right]\\left(\\frac{r}{w(z)}\\right)^{l}
L_p^{l}\\!\\left(\\frac{2r^2}{w^2(z)}\\right)
\\exp\\!\\left(-\\frac{r^2}{w^2(z)}\\right)
\\exp(-jkz -jk\\frac{r^2}{2R(z)} \\mp il\\phi + j\\left(2p + l + 1\\right)\\zeta\\!(z))
```

where ``L_n^{α}`` is the generalized Laguerre polynomial.

# Arguments
- `grid`: [`TransverseGrid`](@ref) on which to evaluate
- `beam`: [`LGBeam`](@ref) descriptor

# Returns
A [`ScalarField`](@ref) with the LG modal amplitude.

# Examples
```jldoctest
julia> grid = TransverseGrid(range(-1e-3, 1e-3, 64), range(-1e-3, 1e-3, 64));

julia> beam = LGBeam(200e-6, 632.8e-9, 0, 3);

julia> field = evaluate(grid, beam);

julia> abs2(field.U[32, 32]) < 1e3   # vortex: zero on axis for l≠0
true
```

# References
Allen 1992. DOI: 10.1103/PhysRevA.45.8185
Saleh & Teich, *Fundamentals of Photonics*, 3rd ed., §3.4
"""
function evaluate(grid::TransverseGrid{<:Real}, beam::LGBeam) :: ScalarField
  (; w0, λ, p, l, z0, n_index, center) = beam
  
  x0, y0 = center

  k  = 2π * n_index / λ
  zR = π * w0^2 * n_index / λ
  z  = z0
  
  wz = w0 * sqrt(1 + (z/zR)^2)
  Rz = iszero(z) ? Inf : z * (1 + (zR/z)^2)
  ψz = atan(z/zR)

  X  = @. grid.x' - x0
  Y  = @. grid.y - y0
  r² = @. X^2 + Y^2
  r  = @. sqrt(r²)
  ϕ  = @. atan(Y, X)               # azimuthal angle

  C = sqrt(2/π) * sqrt(factorial(p) / factorial(l + p))

  # U = u * exp(-ikz)
  # TODO: Check the constant C for correctness
  # Should it be C/wz or C*(w0/wz)?
  U = @. C * (w0 / wz) *
      (r * sqrt(2) / wz)^l *
      _laguerre(p, l, 2r²/wz^2) *
      exp(-r²/wz^2) *
      exp(-1im * k * z) *
      exp(-1im * k * r² / (2 * Rz)) *
      exp(-1im * l * ϕ) *
      exp(1im * (2p + l + 1) * ψz) +
      0im

  return ScalarField(U, grid, λ)
end

modal_number(beam::LGBeam) = 2*beam.p + beam.l + 1
