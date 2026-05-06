"""
    propagate_ssf(field::ScalarField, z::Float64) -> ScalarField

Propagate a scalar field by distance `z` [m] using the Fresnel
single-FFT method. The output grid scales with propagation distance,
making this suitable for diverging beams and large propagation distances.

The output grid spacing is:

```math
\\Delta x_{out} = \\frac{\\lambda z}{N_x \\Delta x_{in}}
```

and the output field is:

```math
E_{out}(x, y) = \\frac{-i}{\\lambda z} \\exp(ikz)
\\exp\\!\\left(\\frac{i\\pi}{\\lambda z}(x^2+y^2)\\right)
\\mathcal{F}\\left\\{E_{in}(x,y)
\\exp\\!\\left(\\frac{i\\pi}{\\lambda z}(x^2+y^2)\\right)\\right\\}
```

# Arguments
- `field`: input `ScalarField`
- `z`    : propagation distance [m]

# Notes
Unlike `propagate`, the output grid is different from the input grid
and grows with `z`. Use this for diverging beams, paraboloidal waves,
or any field that would clip the grid edges with the angular spectrum method.

# References
Goodman, *Introduction to Fourier Optics*, 3rd ed., §4.2
"""
function propagate_ssf(field::ScalarField, z::Float64) :: ScalarField
  # TODO: Single-step Fresnel propagation

  λ = field.λ
  Nx = field.grid.Nx
  Ny = field.grid.Ny
  dx1 = field.grid.dx
  dy1 = field.grid.dy
  Lx1 = Nx*dx1
  Ly1 = Ny*dy1

  @assert dx1 <= λ*z/Lx1 "dx1 should satisfy dx <= λz/Lx"
  @assert dy1 <= λ*z/Ly1 "dx1 should satisfy dx <= λz/Lx"

  k = 2π / λ
  dx2 = λ*z/(Nx*dx1)
  dy2 = λ*z/(Ny*dy1)

  x1 = fftshift(fftfreq(Nx, Nx*dx1))
  y1 = fftshift(fftfreq(Ny, Ny*dy1))
  x2 = fftshift(fftfreq(Nx, Nx*dx2))
  y2 = fftshift(fftfreq(Ny, Ny*dy2))

  Cz = (exp(1im*k*z)/(1im*λ*z))
  quadratic_factor = @. exp(1im*π*(x1'^2 + y1^2)/(λ*z))
  U = @. field.U * quadratic_factor

  U_out = Cz * fftshift(fft(ifftshift(U))) * dx1 * dy1
  grid_out = TransverseGrid(x2, y2)
  λ_out = λ

  return ScalarField(U_out, grid_out, λ_out)
end
