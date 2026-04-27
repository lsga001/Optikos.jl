"""
    mix(crystal::SFGCrystal, field1::ScalarField, field2::ScalarField) -> ScalarField

Apply sum frequency generation mixing between two fields.
Returns a field proportional to E₁ · E₂ with frequency ω₃ = ω₁ + ω₂.

Both fields must share the same grid.
"""
function mix(crystal::SFGCrystal, field1::ScalarField, field2::ScalarField) :: ScalarField
  @assert field1.grid === field2.grid || 
          (field1.grid.Nx == field2.grid.Nx && 
           field1.grid.Ny == field2.grid.Ny &&
           field1.grid.dx ≈ field2.grid.dx &&
           field1.grid.dy ≈ field2.grid.dy) "Fields must share the same grid"

  U_out = @. field1.U * field2.U
  grid_out = field1.grid
  λ_out = (field1.λ^-1 + field2.λ^-1)^-1

  return ScalarField(U_out, grid_out, λ_out)
end
