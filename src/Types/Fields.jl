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
intensity(f::ScalarField) = ScalarField(abs2.(f.U),f.grid,f.λ)
Base.real(f::ScalarField) = ScalarField(complex.(Base.real.(f.U)),f.grid,f.λ)
Base.imag(f::ScalarField) = ScalarField(complex.(imag.(f.U)),f.grid,f.λ)
Base.abs(f::ScalarField) = ScalarField(abs.(f.U),f.grid,f.λ)
Base.sqrt(f::ScalarField) = ScalarField(sqrt.(abs(f).U),f.grid,f.λ)

function Base.:+(f1::ScalarField, f2::ScalarField) :: ScalarField
  @assert f1.λ ≈ f2.λ "Wavelengths must match"
  @assert f1.grid == f2.grid "Grids must match"
  return ScalarField(f1.U .+ f2.U, f1.grid, f1.λ)
end

function Base.:-(f1::ScalarField, f2::ScalarField) :: ScalarField
  @assert f1.λ ≈ f2.λ "Wavelengths must match"
  @assert f1.grid == f2.grid "Grids must match"
  return ScalarField(f1.U .- f2.U, f1.grid, f1.λ)
end

function Base.:*(f1::ScalarField, f2::ScalarField) :: ScalarField
  @assert f1.λ ≈ f2.λ "Wavelengths must match"
  @assert f1.grid == f2.grid "Grids must match"
  return ScalarField(f1.U .* f2.U, f1.grid, f1.λ)
end

function Base.:*(α::Number,f::ScalarField) :: ScalarField
  return ScalarField(α .* f.U, f.grid, f.λ)
end

function Base.:*(f::ScalarField, α::Number) :: ScalarField
  return α * f
end

function Base.:/(f::ScalarField, α::Number) :: ScalarField
  return ScalarField(f.U ./ α, f.grid, f.λ)
end

function Base.conj(f::ScalarField) :: ScalarField
  return ScalarField(conj.(f.U), f.grid, f.λ)
end

#struct 

# TODO: VectorField
# To be added after polarization module is developed

#struct VectorField
#  Ex :: Matrix{ComplexF64}    # Complex Amplitude [V/m]
#  Ey :: Matrix{ComplexF64}    # Complex Amplitude [V/m]
#  grid :: TransverseGrid
#  λ :: Float64                # Wavelength [m]
#end
