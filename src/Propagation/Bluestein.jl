using FFTW

# TODO: Check bluestein propagation implementation.
# We should check that the intensity and amplitudes are correct.

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

  @assert z >= 0 "Propagation distance must be nonnegative"

  λ  = field.λ
  k  = 2π / λ
  xv = field.grid.x'
  yv = field.grid.y
  uv = output_grid.x'
  vv = output_grid.y
  dx = field.grid.dx
  dy = field.grid.dy

  F0 = @. exp(1im*k*z) / (1im*λ*z) * exp(1im*k*(uv^2 + vv^2) / (2z))
  F  = @. exp(1im*π*(xv^2 + yv^2) / (λ*z))

  fxs = λ * z / dx
  fx  = uv .+ fxs / 2
  fys = λ * z / dy
  fy  = vv .+ fys / 2

  E = dx * dy .* F0 .* _fft2_bluestein(field.U .* F, fx, fxs, fy, fys)

  return ScalarField(E, output_grid, λ)
end

function _fft2_bluestein(U, fx, fxs, fy, fys)
    fx1 = fx[1];  fx2 = fx[end];  Mx = length(fx)
    fy1 = fy[1];  fy2 = fy[end];  My = length(fy)
    return _fft_bluestein(
               _fft_bluestein(U,  fx1, fx2, fxs, Mx; axis=1),
               fy1, fy2, fys, My; axis=2
           )
end

function _fft_bluestein(x, f1, f2, fs, M; axis=1)
    x = axis == 1 ? permutedims(x, [1, 2]) : permutedims(x, [2, 1])

    N, n_cols = size(x)

    df      = f2 - f1
    phi1    = 2π * f1 / fs
    phi2    = 2π * f2 / fs
    A       = -exp(1im * phi1)
    W       = exp(-1im * (phi2 - phi1) / M)

    fv      = range(f1, step=df/M, length=M)
    M_shift = -N / 2
    P_shift = @. exp(-1im * 2π * M_shift * fv / fs)

    # Precompute CZT kernel values — same for every column
    k_N   = collect(0:N-1)
    k_M   = collect(0:M-1)
    Wk2_N = @. W^(k_N^2 / 2)    # W^(k²/2) for k=0..N-1
    Wk2_M = @. W^(k_M^2 / 2)    # W^(k²/2) for k=0..M-1
    Ak    = @. A^(-k_N)           # A^(-k)   for k=0..N-1

    # Precompute the Toeplitz structure — same for every column
    n_pad = Int(2^ceil(log2(M + N - 1)))
    r     = conj.(Wk2_N)
    c     = conj.(Wk2_M)
    ft    = _precompute_toeplitz(r, c, N, M, n_pad)

    # Apply CZT to each column reusing precomputed values
    X = zeros(ComplexF64, M, n_cols)
    for j in 1:n_cols
        X[:, j] = P_shift .* _czt(x[:, j], M, Wk2_N, Wk2_M, Ak, ft, n_pad)
    end

    return axis == 1 ? permutedims(X, [1, 2]) : permutedims(X, [2, 1])
end


"""
    _precompute_toeplitz(r, c, N, M, n_pad) -> Vector{ComplexF64}

Precompute the FFT of the Toeplitz matrix embedding, which is the same
for every column and can be reused across all CZT calls with the same
parameters.
"""
function _precompute_toeplitz(r, c, N, M, n_pad)
    c_hat = zeros(ComplexF64, n_pad)
    c_hat[1:M] .= c
    for i in 2:N
        c_hat[n_pad - i + 2] = r[i]
    end
    return fft(c_hat)
end


"""
    _czt(x, M, Wk2_N, Wk2_M, Ak, ft, n_pad) -> Vector{ComplexF64}

Compute the chirp-z transform of vector `x` using precomputed kernel
values. Accepts precomputed `Wk2_N`, `Wk2_M`, `Ak`, and the FFT of the
Toeplitz embedding `ft` to avoid redundant computation across columns.
"""
function _czt(x, M, Wk2_N, Wk2_M, Ak, ft, n_pad)
    N = length(x)

    # Step 1: multiply input by Bluestein kernel
    X = @. Wk2_N * Ak * x

    # Step 2: Toeplitz multiplication via circular convolution
    X = _toeplitz_multiply(X, ft, N, M, n_pad)

    # Step 3: multiply output by Bluestein kernel
    X = @. Wk2_M * X

    return X
end


"""
    _toeplitz_multiply(x, ft, N, M, n_pad) -> Vector{ComplexF64}

Compute the Toeplitz matrix-vector product via circular convolution,
using a precomputed FFT of the Toeplitz embedding `ft`.
"""
function _toeplitz_multiply(x, ft, N, M, n_pad)
    x_hat        = zeros(ComplexF64, n_pad)
    x_hat[1:N]  .= x
    y_hat        = _circulant_multiply(ft, x_hat)
    return y_hat[1:M]
end


"""
    _circulant_multiply(ft, x) -> Vector{ComplexF64}

Compute the circulant matrix-vector product y = Gx using the
precomputed FFT of the circulant generator `ft`.
"""
function _circulant_multiply(ft, x)
    return ifft(ft .* fft(x))
end

#function _fft2_bluestein(U, fx, fxs, fy, fys)
#  fx1 = fx[1];   fx2 = fx[end];   Mx = length(fx)
#  fy1 = fy[1];   fy2 = fy[end];   My = length(fy)
#  return _fft_bluestein(
#           _fft_bluestein(U, fx1, fx2, fxs, Mx; axis=1),
#           fy1, fy2, fys, My; axis=2
#         )
#end
#
#function _fft_bluestein(x, f1, f2, fs, M; axis=1)
#  x = axis == 1 ? permutedims(x, [1, 2]) : permutedims(x, [2, 1])
#
#  N, n_cols = size(x)
#
#  df = f2-f1
#
#  phi1 = 2π * f1 / fs
#  phi2 = 2π * f2 / fs
#
#  A = -exp(1im * phi1)
#  W = exp(-1im * (phi2 - phi1) / M)
#
#  M_shift = -N/2   
#  fv  = range(f1, step=df/M, length=M)              # [f1, f1 + (1/M)df, ..., f1 + ((M-1)/M) * df]
#  P_shift = @. exp(-1im * 2π * M_shift * fv / fs)
#
#  X = mapslices(u -> P_shift .* _czt(u, M, A, W), x; dims=1)
#
#  X = axis == 1 ? permutedims(X, [1, 2]) : permutedims(X, [2, 1])
#
#  return X
#end
#
#
#function _czt(x, M, A, W)
#  N = size(x, 1)
#  X = zeros(ComplexF64, N)
#  r = zeros(ComplexF64, N)
#  c = zeros(ComplexF64, M)
#
#  for i in 1:N
#    k    = i - 1
#    X[i] = W^(k^2/2) * A^(-k) * x[i]
#    r[i] = W^(-k^2/2)
#  end
#
#  for i in 1:M
#    k    = i - 1
#    c[i] = W^(-k^2/2)
#  end
#
#  X = _toeplitz_multiply(r, c, X)
#
#  for i in 1:M
#    k    = i - 1
#    X[i] = W^(k^2/2) * X[i]
#  end
#
#  return X
#end
#
#
#function _toeplitz_multiply(r, c, x)
#  N = length(r)
#  M = length(c)
#  n = Int(2^ceil(log2(M + N - 1)))
#
#  c_hat = zeros(ComplexF64, n)
#  for i in 1:M
#    c_hat[i] = c[i]
#  end
#  for i in 2:N
#    c_hat[n - i + 2] = r[i]
#  end
#
#  x_hat = zeros(ComplexF64, n)
#  for i in 1:N
#    x_hat[i] = x[i]
#  end
#
#  y_hat = _circulant_multiply(c_hat, x_hat)
#
#  y = zeros(ComplexF64, M)
#  for i in 1:M
#    y[i] = y_hat[i]
#  end
#
#  return y
#end
#
#
#function _circulant_multiply(c, x)
#  C = fft(c)
#  X = fft(x)
#  Y = zeros(ComplexF64, length(c))
#  for i in 1:length(c)
#    Y[i] = C[i] * X[i]
#  end
#  return ifft(Y)
#end
#
## TODO: Need to figure out where to use zeropad
#function _zeropad(x, n)
#  m     = length(x)
#  x_hat = zeros(ComplexF64, n)
#  for i in 1:m
#    x_hat[i] = x[i]
#  end
#  return x_hat
#end
