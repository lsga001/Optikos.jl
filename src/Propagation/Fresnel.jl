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
  @assert dy1 <= λ*z/Ly1 "dy1 should satisfy dy <= λz/Ly"

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


"""
"Algorithm for reconstruction of digital holograms with adjustable magnification" - 
Fucai Zhang, Ichirou Yamaguchi, and L. P. Yaroslavsky
DOI: 10.1364/OL.29.001668
"""
function propagate_dsf(field::ScalarField, z::Float64, m::Real) :: ScalarField
  # TODO: Double-step Fresnel propagation

  λ = field.λ
  Nx = field.grid.Nx
  Ny = field.grid.Ny
  dx1 = field.grid.dx
  dy1 = field.grid.dy
  Lx1 = Nx*dx1
  Ly1 = Ny*dy1

  k = 2π / λ
  
  # This relationship is satisfied: z = z1 + z2
  
  z1 = z / (1 + m)
  z2 = z * (+m) / (1 + m)

  ft1 = (z1 >= 0) ? fft : ifft
  ft2 = (z2 >= 0) ? fft : ifft
  
  dxv = λ*z1 / (Nx * dx1)
  dyv = λ*z1 / (Ny * dy1)

  dx2 = λ*z2 / (Nx * dxv)
  dy2 = λ*z2 / (Ny * dyv)

  x1 = fftshift(fftfreq(Nx, Nx*dx1))
  y1 = fftshift(fftfreq(Ny, Ny*dy1))
  
  xv = fftshift(fftfreq(Nx, Nx*dxv))
  yv = fftshift(fftfreq(Ny, Ny*dyv))

  x2 = fftshift(fftfreq(Nx, Nx*dx2))
  y2 = fftshift(fftfreq(Ny, Ny*dy2))
  
  @assert dx1 <= abs(λ*z1/Lx1) "dx1 should satisfy dx1 <= λz1/Lx1"
  @assert dy1 <= abs(λ*z1/Ly1) "dy1 should satisfy dy1 <= λz1/Ly1"

  xv_max = λ*z1*z2 / (2*z*dxv)
  yv_max = λ*z1*z2 / (2*z*dyv)

  Cz1 = (exp(1im*k*z1)/(1im*λ*z1))
  Cz2 = (exp(1im*k*z2)/(1im*λ*z2))
  Cz = (exp(1im*k*z)/(1im*λ*z))

  quadratic_factor = @. exp(1im*π*(x1'^2 + y1^2)/(λ*z1))
  U = @. field.U * quadratic_factor

  U_1 = Cz1 * fftshift(ft1(ifftshift(U))) * dx1 * dy1

  chirp = @. exp(1im*π*z*(xv'^2 + yv^2)/(λ*z1*z2))
  rect = @. _rect(xv'/xv_max) * _rect(yv/yv_max)

  U_2 = fftshift(ft2(ifftshift( chirp .* rect .* U_1 ))) * dxv * dyv

  U_out = Cz2 * U_2

  grid_out = TransverseGrid(x2, y2)
  λ_out = λ

  return ScalarField(U_out, grid_out, λ_out)
end


"""
Far-field optical propagators with user-defined object-plane pixel size for ptychography.
ANTONIOS PELEKANIDIS,1,2 KJELD S. E. EIKEMA,1,2 AND STEFAN WITTE
DOI: 10.1364/OPTCON.555163
"""
function propagate_2sf(field::ScalarField, z::Float64, x_max::Float64)
  # TODO: Unfinished two step fresnel rework

  λ = field.λ
  Nx = field.grid.Nx
  Ny = field.grid.Ny
  dx1 = field.grid.dx
  dy1 = field.grid.dy
  Lx1 = Nx*dx1
  Ly1 = Ny*dy1
  Lx2 = 2 * x_max

  dq = 2*x_max/Nx



end

function _rect(t::Real)
  t_abs = abs(t)
  Result =  (t_abs > 1/2) ? 0 :
            (t_abs == 1/2) ? 1/2 : 1
  return Result
end
