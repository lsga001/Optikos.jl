using FFTW

# TODO: Check bluestein propagation implementation.
"""
    propagate_bluestein(field::ScalarField, z::Float64,
                        output_grid::TransverseGrid) -> ScalarField

Scalar diffraction propagation using the Bluestein/chirp-z transform method.
Allows independent control of input and output grid size and spacing,
avoiding the wraparound artifacts of standard FFT-based propagators.

The algorithm computes:

```math
E_{out}(u, v) = \\frac{e^{ikz}}{i\\lambda z}
\\exp\\!\\left(\\frac{ik(u^2+v^2)}{2z}\\right)
\\Delta x \\Delta y \\,
\\mathcal{F}_{\\text{CZT}}\\left\\{E_{in}(x,y)
\\exp\\!\\left(\\frac{i\\pi(x^2+y^2)}{\\lambda z}\\right)\\right\\}
```

where ``\\mathcal{F}_{\\text{CZT}}`` is the chirp-z transform evaluated
on the output grid coordinates.

# Arguments
- `field`      : input `ScalarField` at the source plane
- `z`          : propagation distance [m]
- `output_grid`: `TransverseGrid` defining the observation plane sampling.
                 Can have different size and spacing from the input grid.

# Notes
This method is particularly suited for:
- Diverging beams (paraboloidal, spherical waves)
- Large propagation distances where the field spreads beyond the input grid
- Cases requiring different input/output resolutions

# References
Hu et al., *Light: Science and Applications* 9, 119 (2020).
DOI: 10.1038/s41377-020-00362-z
"""
function propagate_bluestein(
    field       :: ScalarField,
    z           :: Float64,
    output_grid :: TransverseGrid{<:Real},
) :: ScalarField

    @assert z > 0 "Propagation distance must be positive"

    λ  = field.λ
    k  = 2π / λ

    xv = field.grid.x
    yv = field.grid.y
    uv = output_grid.x
    vv = output_grid.y

    dx = field.grid.dx
    dy = field.grid.dy

    F0 = @. exp(1im*k*z) / (1im*λ*z) * exp(1im*k*(uv^2 + vv'^2) / (2z))
    F  = @. exp(1im*π*(xv^2 + yv'^2) / (λ*z))

    fxs = λ * z / dx
    fx  = uv .+ fxs / 2
    fys = λ * z / dy
    fy  = vv .+ fys / 2

    E = dx * dy .* F0 .* _fft2_bluestein(field.U .* F, fx, fxs, fy, fys)

    return ScalarField(E, output_grid, λ)
end


function _fft2_bluestein(U, fx, fxs, fy, fys)
    fx1 = fx[1];   fx2 = fx[end];   Mx = length(fx)
    fy1 = fy[1];   fy2 = fy[end];   My = length(fy)
    return _fft_bluestein(
               _fft_bluestein(U, fx1, fx2, fxs, Mx; axis=1),
               fy1, fy2, fys, My; axis=2
           )
end


function _fft_bluestein(x, f1, f2, fs, M; axis=1)
    if axis == 1
        x = permutedims(x, [1, 2])
    else
        x = permutedims(x, [2, 1])
    end

    N = size(x, 1)

    f11 = f1 + (M*fs + f2 - f1) / (2M)
    f22 = f2 + (M*fs + f2 - f1) / (2M)

    phi1 = 2π * f11 / fs
    phi2 = 2π * f22 / fs

    A = exp(1im * phi1)
    W = exp(-1im * (phi2 - phi1) / M)

    fv      = range(f11, f22, M)
    M_shift = -N/2 + 1/2
    P_shift = @. exp(-1im * 2π * fv * M_shift / fs)

    X = mapslices(u -> P_shift .* _czt(u, M, A, W), x; dims=1)

    if axis == 1
        X = permutedims(X, [1, 2])
    else
        X = permutedims(X, [2, 1])
    end

    return X
end


function _czt(x, M, A, W)
    N = size(x, 1)
    X = zeros(ComplexF64, N)
    r = zeros(ComplexF64, N)
    c = zeros(ComplexF64, M)

    for i in 1:N
        k    = i - 1
        X[i] = W^(k^2/2) * A^(-k) * x[i]
        r[i] = W^(-k^2/2)
    end

    for i in 1:M
        k    = i - 1
        c[i] = W^(-k^2/2)
    end

    X = _toeplitz_multiply(r, c, X)

    for i in 1:M
        k    = i - 1
        X[i] = W^(k^2/2) * X[i]
    end

    return X
end


function _toeplitz_multiply(r, c, x)
    N = length(r)
    M = length(c)
    n = Int(2^ceil(log2(M + N - 1)))

    c_hat = zeros(ComplexF64, n)
    for i in 1:M
        c_hat[i] = c[i]
    end
    for i in 2:N
        c_hat[n - i + 2] = r[i]
    end

    x_hat = zeros(ComplexF64, n)
    for i in 1:N
        x_hat[i] = x[i]
    end

    y_hat = _circulant_multiply(c_hat, x_hat)

    y = zeros(ComplexF64, M)
    for i in 1:M
        y[i] = y_hat[i]
    end

    return y
end


function _circulant_multiply(c, x)
    C = fft(c)
    X = fft(x)
    Y = zeros(ComplexF64, length(c))
    for i in 1:length(c)
        Y[i] = C[i] * X[i]
    end
    return ifft(Y)
end


function _zeropad(x, n)
    m     = length(x)
    x_hat = zeros(ComplexF64, n)
    for i in 1:m
        x_hat[i] = x[i]
    end
    return x_hat
end
