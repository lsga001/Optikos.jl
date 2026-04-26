"""
    propagate_fresnel(field::ScalarField, z::Float64) -> ScalarField

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
function propagate_fresnel(field::ScalarField, z::Float64) :: ScalarField
    @assert z > 0 "Propagation distance must be positive"

    k  = 2π / field.λ
    Nx = field.grid.Nx
    Ny = field.grid.Ny
    dx = field.grid.dx
    dy = field.grid.dy

    # Output grid spacing — scales with z
    dx_out = field.λ * z / (Nx * dx)
    dy_out = field.λ * z / (Ny * dy)

    # Output grid coordinates
    x_out = range(-(Nx÷2) * dx_out, (Nx÷2 - 1) * dx_out, length=Nx)
    y_out = range(-(Ny÷2) * dy_out, (Ny÷2 - 1) * dy_out, length=Ny)

    # Input coordinates for quadratic phase
    X_in = field.grid.x'
    Y_in = field.grid.y

    # Output coordinates for quadratic phase
    X_out = x_out'
    Y_out = y_out

    # Fresnel propagator — single FFT method
    # 1. Multiply input by quadratic phase
    quadratic_in = @. exp(1im * k / (2z) * (X_in^2 + Y_in^2))

    # 2. FFT
    E_fft = fftshift(fft(ifftshift(field.U .* quadratic_in)))

    # 3. Multiply by output quadratic phase and prefactor
    prefactor = -1im / (field.λ * z) * exp(1im * k * z) * dx * dy
    quadratic_out = @. exp(1im * k / (2z) * (X_out^2 + Y_out^2))

    E_out = prefactor .* quadratic_out .* E_fft

    new_grid = TransverseGrid(x_out, y_out)
    return ScalarField(E_out, new_grid, field.λ)
end
