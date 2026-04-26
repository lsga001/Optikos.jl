"""
    apply(space::FreeSpace, field::ScalarField) -> ScalarField

Propagates a field a distance z defined by `space::FreeSpace`.
"""
function apply(space::FreeSpace, field::ScalarField) :: ScalarField
  return propagate_fresnel(field, space.z)
end

function apply(lens::ThinLens, field::ScalarField) :: ScalarField
  k  = 2π / field.λ
  X  = field.grid.x'
  Y  = field.grid.y
  r² = @. X^2 + Y^2
  t  = @. exp(-1im * k * r² / (2 * lens.f))
  return ScalarField(field.U .* t, field.grid, field.λ)
end

function apply(aperture::CircularAperture, field::ScalarField) :: ScalarField
  X = field.grid.x'
  Y = field.grid.y
  mask = @. (X^2 + Y^2 <= aperture.R^2)
  return ScalarField(field.U .* mask, field.grid, field.λ)
end

function apply(mask::PhaseMask, field::ScalarField) :: ScalarField
  @assert size(mask.φ) == size(field.U) "Phase mask size must match field size"
  t = @. exp(1im * mask.φ)
  return ScalarField(field.U .* t, field.grid, field.λ)
end
