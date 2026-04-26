"""
    DiskSchellModel(R, λ)

Analytical descriptor of a spatially partially coherent source with a
uniform circular emission profile of radius `R`.

Each realization is a paraboloidal wave emitted from a point drawn
uniformly at random from within a disk of radius `R`. The degree of
spatial coherence is controlled by `R`:

- `R → 0` : approaches a fully coherent point source
- `R → ∞` : approaches a fully spatially incoherent extended source

The coherence width of the resulting field is approximately `λz/R`
after propagating a distance `z`, by the van Cittert-Zernike theorem.

# Fields
- `R`: source disk radius [m]
- `λ`: wavelength [m]

# Examples
```jldoctest
julia> source = DiskSchellModel(500e-6, 632.8e-9);

julia> source.R
0.0005
```

# References
Mandel & Wolf, *Optical Coherence and Quantum Optics*, §4.4
Goodman, *Statistical Optics*, 2nd ed., §5.4 (van Cittert-Zernike theorem)
"""
struct DiskSchellModel <: AbstractEnsembleSource
    R :: Float64
    λ :: Float64

    function DiskSchellModel(R, λ; n=1.0)
        @assert R > 0 "Source radius must be positive"
        @assert λ > 0 "Wavelength must be positive"
        new(R, λ)
    end
end

wavelength(s::DiskSchellModel) = s.λ

# TODO: Check the behavior of this type
"""
    GaussianSchellModel(σ_I, σ_g, λ; A=1.0)

Analytical descriptor of a Gaussian Schell-model (GSM) partially coherent source.

The cross-spectral density takes the form:

```math
W(\\mathbf{r}_1, \\mathbf{r}_2) = A^2
\\exp\\!\\left(-\\frac{r_1^2 + r_2^2}{4\\sigma_I^2}\\right)
\\exp\\!\\left(-\\frac{|\\mathbf{r}_1 - \\mathbf{r}_2|^2}{2\\sigma_g^2}\\right)
```

Each realization is a Gaussian beam with unit waist `σ_I`, displaced by
a random amount drawn from a Gaussian distribution of width `σ_g`. The
coherence width `σ_g` controls the degree of spatial coherence:

- `σ_g → ∞` : approaches a fully coherent Gaussian beam
- `σ_g → 0` : approaches a spatially incoherent source

# Fields
- `A`  : field amplitude [V/m]
- `σ_I`: intensity width (rms beam radius) [m]
- `σ_g`: coherence width [m]
- `λ`  : wavelength [m]

# Examples
```jldoctest
julia> source = GaussianSchellModel(500e-6, 100e-6, 632.8e-9);

julia> source.σ_g < source.σ_I    # partially coherent
true
```

# References
Mandel & Wolf, *Optical Coherence and Quantum Optics*, §5.6
Goodman, *Statistical Optics*, 2nd ed., §5.6
"""
struct GaussianSchellModel <: AbstractEnsembleSource
    A   :: Float64
    σ_I :: Float64
    σ_g :: Float64
    λ   :: Float64

    function GaussianSchellModel(σ_I, σ_g, λ; A=1.0)
        @assert σ_I > 0 "Intensity width must be positive"
        @assert σ_g > 0 "Coherence width must be positive"
        @assert λ   > 0 "Wavelength must be positive"
        new(A, σ_I, σ_g, λ)
    end
end

wavelength(s::GaussianSchellModel) = s.λ
