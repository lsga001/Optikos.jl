"""
    ScalarField(U::Matrix{ComplexF64}, grid::TransverseGrid, λ::Float64)

A transverse scalar field representation type.

# Arguments
- `U`: Complex amplitude [V/m]
- `grid`: A uniform transverse grid associated with the complex amplitude
- `λ`: Wavelength [m]

# Examples

```jldoctest
julia> x = y = range(-1e-3, 1e-3, length=256);

julia> grid = TransverseGrid(x, y);

julia> λ = 600e-9;

julia> U = ones(ComplexF64, (256, 256));

julia> field = ScalarField(U, grid, λ);
```
"""
struct ScalarField <: AbstractField
  U :: Matrix{ComplexF64}     # Complex Amplitude [V/m]
  grid :: TransverseGrid{<:Real}
  λ :: Float64                # Wavelength [m]

  function ScalarField(U, grid::TransverseGrid{<:Real}, λ)
    @assert size(U) == (grid.Nx, grid.Ny) "Field matrix size must match grid dimensions"
    @assert λ > 0 "Wavelength must be positive"
    new(U, grid, λ)
  end
end

# Fulfill AbstractField interface
grid(f::ScalarField) = f.grid
wavelength(f::ScalarField) = f.λ
intensity(f::ScalarField) = Base.abs2.(f.U)
real(f::ScalarField) = Base.real.(f.U)
imaginary(f::ScalarField) = Base.imag.(f.U)
abs(f::ScalarField) = Base.abs.(f.U)



#struct 

# TODO: VectorField
# To be added after polarization module is developed

#struct VectorField
#  Ex :: Matrix{ComplexF64}    # Complex Amplitude [V/m]
#  Ey :: Matrix{ComplexF64}    # Complex Amplitude [V/m]
#  grid :: TransverseGrid
#  λ :: Float64                # Wavelength [m]
#end
