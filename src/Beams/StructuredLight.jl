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
    evaluate(grid::TransverseGrid, beam::GaussianBeam) -> ScalarField

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
function evaluate(grid::TransverseGrid, beam::GaussianBeam) :: ScalarField
  (; w0, λ, z0, n_index) = beam

  k = 2π * n_index / λ
  zR = π * w0^2 * n_index / λ # Rayleigh range
  z = z0
  
  wz = w0 * sqrt(1 + (z/zR)^2 ) # beam waist at z
  Rz = iszero(z) ? Inf : z * (1 + (zR / z)^2) # wavefront radius at z
  ψz = atan(z / zR) # Gouy phase

  X = grid.x'
  Y = grid.y
  r² = @. X^2 + Y^2

  C = sqrt(2/π) / wz

  U = @. C *
      exp(-r² / wz^2) *
      exp(-1im * (k*z + k*r² / (2*Rz) - ψz))

  return ScalarField(U, grid, λ)
end


"""
    evaluate(grid::TransverseGrid, beam::LGBeam) -> ScalarField

Sample a Laguerre-Gaussian beam LG(p, l) onto a transverse grid at `z = z0`.

The complex field amplitude at the waist is:

```math
U(r, \\phi) = \\left(\\frac{r\\sqrt{2}}{w_0}\\right)^{|l|}
L_p^{|l|}\\!\\left(\\frac{2r^2}{w_0^2}\\right)
\\exp\\!\\left(-\\frac{r^2}{w_0^2}\\right)
\\exp(-il\\phi)
```

where ``L_p^{|l|}`` is the generalized Laguerre polynomial.

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
Allen et al., *Phys. Rev. A* 45, 8185 (1992)
Saleh & Teich, *Fundamentals of Photonics*, 3rd ed., §3.3
"""
function evaluate(grid::TransverseGrid, beam::LGBeam) :: ScalarField
    (; w0, λ, p, l, z0, n_index) = beam

    k  = 2π * n_index / λ
    zR = π * w0^2 * n_index / λ
    z  = z0
    abs_l = abs(l)
    
    wz = w0 * sqrt(1 + (z/zR)^2)
    Rz = iszero(z) ? Inf : z * (1 + (zR/z)^2)
    ψz = (2p + abs_l + 1)*atan(z/zR)

    X  = grid.x'
    Y  = grid.y
    r² = @. X^2 + Y^2
    r  = @. sqrt(r²)
    ϕ  = @. atan(Y, X)               # azimuthal angle

    C = sqrt(2 * factorial(p) / (π * factorial(abs_l + p))) / wz
    ρ  = @. r * sqrt(2) / wz

    U = @. C *
        (ρ^abs_l) *
        _laguerre(p, abs_l, 2r²/wz^2) *
        exp(-r²/wz^2) *
        exp(-1im * k * z) *
        exp(-1im * k * r² / (2 * Rz)) *
        exp(1im * ψz) *
        exp(-1im * l * ϕ) + 0im

    return ScalarField(U, grid, λ)
end

modal_number(beam::LGBeam) = 2*beam.p + abs(beam.l) + 1
