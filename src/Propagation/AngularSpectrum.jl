# TODO: CHECK accuracy
"""
    propagate_angular(field::ScalarField, z::Float64) -> ScalarField

Propagate a scalar field by distance `z` [m] using the angular spectrum method.

The transfer function in the paraxial approximation is:

```math
H(k_x, k_y, z) = \\exp(ikz) \\exp\\!\\left(-iz\\frac{k_x^2 + k_y^2}{2k}\\right)
```

The field is transformed to the spatial frequency domain via FFT,
multiplied by the transfer function, and transformed back.

# Arguments
- `field`: input `ScalarField` at z=0
- `z`    : propagation distance [m], must be positive

# Notes
Assumes uniform grid spacing. The grid is preserved — the output
field lives on the same transverse grid as the input.

# Examples
```jldoctest
julia> grid  = TransverseGrid(range(-2e-3, 2e-3, 256), range(-2e-3, 2e-3, 256));

julia> beam  = GaussianBeam(200e-6, 632.8e-9);

julia> field = evaluate(grid, beam);

julia> field_z = propagate_angular(field, 0.1);
```

# References
Goodman, *Introduction to Fourier Optics*, 3rd ed., §3.10
Saleh & Teich, *Fundamentals of Photonics*, 3rd ed., §4.4
"""
function propagate_angular(field::ScalarField, z::Float64) :: ScalarField
  @assert z >= 0 "Propagation distance must be nonnegative"

  Nx = field.grid.Nx
  Ny = field.grid.Ny
  dx = field.grid.dx
  dy = field.grid.dy
  Lx = Nx*dx/2
  Ly = Nx*dy/2

  fx = fftshift(fftfreq(Nx, 1/dx))'   # shape (1, Nx)
  fy = fftshift(fftfreq(Ny, 1/dy))    # shape (Ny, 1)
  #fx = collect(-1/(2dx):1/(2Lx):1/(2dx)-1/(2Lx))'
  #fy = collect(-1/(2dy):1/(2Ly):1/(2dy)-1/(2Ly))

  λ = field.λ
  k = 2π/λ

  kz² = @. (2π)^2 * (λ^-2 - fx^2 - fy^2)
  kz = @. sqrt(kz²)
  #kz[kz² .>= 0] = kz[kz² .>= 0]
  #kz[kz² .< 0] = -1im*kz[kz² .< 0]

  # Transfer function (Angular spectrum)
  H = @. exp(-1im * kz * z)
  
  # Transfer function (Fresnel approximation)
  #H = @. exp(-1im * k * z) * exp(1im*π*λ*z*(fx^2 + fy^2))
 
  # Angular spectrum propagation
  U = fft(field.U)
  
  U = fftshift(U)

  E_propagated = ifft(ifftshift(U .* H))

  return ScalarField(E_propagated, field.grid, field.λ)
end
