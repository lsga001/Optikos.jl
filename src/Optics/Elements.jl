"""
    apply(space::FreeSpace, field::ScalarField) -> ScalarField

Propagates a field a distance z defined by `space::FreeSpace`.
"""
function apply(space::FreeSpace, field::ScalarField) :: ScalarField
  return propagate_angular(field, space.z)
  #return propagate_fresnel(field, space.z)
  #grid_out = (space.grid_out !== missing) ? space.grid_out : field.grid
  #return propagate_bluestein(field, space.z, grid_out)
end

function apply(lens::ThinLens, field::ScalarField) :: ScalarField 
  f = lens.f
  D_L = lens.diameter
  d0 = lens.thickness
  f_number = abs(f)/D_L

  λ = field.λ
  dx = field.grid.dx
  dy = field.grid.dy

  @assert f_number >= dx/λ "f-number must be bigger than Δx/λ"
  @assert f_number >= dy/λ "f-number must be bigger than Δy/λ"

  k  = 2π / field.λ
  X  = field.grid.x'
  Y  = field.grid.y
  R² = @. X^2 + Y^2
  h₀ = @. exp(im*k*d0) # constant prefactor due to lens thickness

  P = @. (R² <= (D_L/2)^2) # lens pupil function
  t  = @. h₀ * P * exp(-1im * k * R² / (2 * lens.f))
  return ScalarField(field.U .* t, field.grid, field.λ)
end

function apply(lens::FourierLens, field::ScalarField) :: ScalarField
  # Assumes waves are paraxial and the Fresnel approximation
  λ = field.λ
  k = 2π/λ
  d = lens.d # Distance before lens
  f = lens.f # Focal length of lens

  X = field.grid.x'
  Y = field.grid.y

  dx = field.grid.dx
  dy = field.grid.dy
  Nx = field.grid.Nx
  Ny = field.grid.Ny
 
  fx = fftshift(fftfreq(Nx, 1 / dx))
  fy = fftshift(fftfreq(Ny, 1 / dy))

  x_out = (λ*f) .* fx
  y_out = (λ*f) .* fy

  hl = (1im/(λ*f)) * exp(-1im*k*(d+f)) * dx * dy
  prefactor = @. hl * exp(1im*π*(X^2 + Y^2)*(d-f)/(λ*f^2))

  U_out = prefactor .* fftshift(fft(ifftshift(field.U)))
  grid_out = TransverseGrid(x_out, y_out)
  λ_out = λ

  return ScalarField(U_out, grid_out, λ_out)
end

function apply(aperture::CircularAperture, field::ScalarField) :: ScalarField
  X = field.grid.x'
  Y = field.grid.y
  mask = @. (X^2 + Y^2 <= aperture.R^2)
  return ScalarField(field.U .* mask, field.grid, field.λ)
end

function apply(spe::SpiralPhaseElement, field::ScalarField) :: ScalarField
  X = field.grid.x'
  Y = field.grid.y
  φ = @. atan(Y, X)
  t = @. cis(spe.l * φ)
  return ScalarField(field.U .* t, field.grid, field.λ)
end

function apply(::Rotate180, field::ScalarField) :: ScalarField
  #@assert field.grid.x[1] ≈ -field.grid.x[end] "Grid must be symmetric around zero"
  #@assert field.grid.y[1] ≈ -field.grid.y[end] "Grid must be symmetric around zero"

  U_out = rot180(field.U) # Rotate by 180. U[i,j] -> U[M-i, N-j]
  grid_out = field.grid
  λ_out = field.λ

  return ScalarField(U_out, grid_out, λ_out)
end

function apply(::Conjugate, field::ScalarField) :: ScalarField
  U_out = @. conj(field.U) # Rotate by 180. U[i,j] -> U[M-i, N-j]
  grid_out = field.grid
  λ_out = field.λ

  return ScalarField(U_out, grid_out, λ_out)
end

function apply(::ConjugateInverter, field::ScalarField) :: ScalarField
  #@assert field.grid.x[1] ≈ -field.grid.x[end] "Grid must be symmetric around zero"
  #@assert field.grid.y[1] ≈ -field.grid.y[end] "Grid must be symmetric around zero"

  U_out = conj.(rot180(field.U)) # Rotate by 180. U[i,j] -> U*[M-i, N-j]
  grid_out = field.grid
  λ_out = field.λ

  return ScalarField(U_out, grid_out, λ_out)
end

function apply(slm::SLM, field::ScalarField) :: ScalarField
  k = 2π/field.λ
  X = field.grid.x'
  Y = field.grid.y

  kx, ky, = get_wavevectors(k, slm.θ_xy, slm.θ_z)

  t = @. exp(1im * slm.φ + kx*X + ky*Y)

  U_out = @. field.U * t
  grid_out = field.grid
  λ_out = field.λ

  return ScalarField(U_out, grid_out, λ_out)
end

function apply(mask::PhaseMask, field::ScalarField) :: ScalarField
  @assert size(mask.φ) == size(field.U) "Phase mask size must match field size"
  t = @. exp(1im * mask.φ)
  return ScalarField(field.U .* t, field.grid, field.λ)
end

