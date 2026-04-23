struct ScalarField <: AbstractField
  E :: Matrix{ComplexF64}     # Complex Amplitude [V/m]
  λ :: Float64                # Wavelength [m]
  grid :: TransverseGrid
end

grid(f::ScalarField) = f.grid
wavelength(f::ScalarField) = f.λ
intensity(f::ScalarField) = abs2.(f.E)

#struct 

# TODO: VectorField
# To be added after polarization module is developed

#struct VectorField
#  Ex :: Matrix{ComplexF64}    # Complex Amplitude [V/m]
#  Ey :: Matrix{ComplexF64}    # Complex Amplitude [V/m]
#  grid :: TransverseGrid
#  λ :: Float64                # Wavelength [m]
#end
