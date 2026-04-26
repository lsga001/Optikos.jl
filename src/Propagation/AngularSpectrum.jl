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

    k  = 2π / field.λ
    Nx = field.grid.Nx
    Ny = field.grid.Ny
    dx = field.grid.dx
    dy = field.grid.dy

    # Spatial frequency grids [rad/m]
    # fftfreq returns frequencies in cycles/sample, multiply by 2π/dx for rad/m
    kx = fftfreq(Nx, 2π / dx)'   # shape (1, Nx)
    ky = fftfreq(Ny, 2π / dy)    # shape (Ny, 1)

    # Paraxial transfer function
    H = @. exp(-1im * z * (kx^2 + ky^2) / (2k)) * exp(1im * k * z)

    # Angular spectrum propagation
    E_propagated = ifft(fft(field.U) .* H)

    return ScalarField(E_propagated, field.grid, field.λ)
end
