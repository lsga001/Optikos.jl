"""
    FreeSpace(z)

Free space propagation by distance `z` [m].
"""
struct FreeSpace <: AbstractOpticalElement
  z :: Float64
  grid_out :: Union{Missing, TransverseGrid{<:Real}}
end

function FreeSpace(z; grid_out=missing)
  @assert z >= 0 "Propagation distance must be nonnegative"
  return FreeSpace(z, grid_out)
end

"""
    ThinLens(f)

Thin lens with focal length `f` [m].
Applies a quadratic phase mask:

```math
t(x,y) = \\exp\\!\\left(-i\\frac{k}{2f}(x^2 + y^2)\\right)
```

Positive `f` for converging, negative for diverging.
"""
struct ThinLens <: AbstractOpticalElement
    f :: Float64

    function ThinLens(f::Float64)
        #@assert !iszero(f) "Focal length cannot be zero"
        new(f)
    end
end

"""
    FourierLens(f)

Fourier lens with focal length `f` [m].
Takes the field to the far-field using a lens with 
focal length `f`.
"""
struct FourierLens <: AbstractOpticalElement
    d :: Float64
    f :: Float64

    function FourierLens(d::Float64, f::Float64)
        @assert d>=0 "Distance before lens cannot be negative"
        @assert !iszero(f) "Focal length cannot be zero"
        new(d,f)
    end
end


"""
    CircularAperture(R)

Hard circular aperture of radius `R>=0` [m].
Transmits the field within radius R and blocks everything outside.
"""
struct CircularAperture <: AbstractOpticalElement
    R :: Float64

    function CircularAperture(R::Float64)
        @assert R >= 0 "Aperture radius must be nonnegative"
        new(R)
    end
end

"""
    SpiralPhaseElement(l)

A spiral phase element (SPE) that imprints an azimuthal phase of
`exp(ilφ)` on an incident field, where `l` is the topological charge
and `φ = atan(y, x)` is the azimuthal angle.

Used to convert a Gaussian or spherical wave into a vortex beam
carrying orbital angular momentum of `lħ` per photon.

# Fields
- `l`: topological charge (integer, positive or negative)

# Examples
```jldoctest
julia> spe = SpiralPhaseElement(1);

julia> spe.l
1
```

# References
Allen et al., *Phys. Rev. A* 45, 8185 (1992)
"""
struct SpiralPhaseElement <: AbstractOpticalElement
    l :: Int
end

"""
    Rotate180()

Optical element that computes 180 degree rotation of the
field: E(x,y) → E(-x,-y).

Physically, this is implemented by a beam inverter (Dove prism
or similar).

Requires the field grid to be symmetric around zero.
"""
struct Rotate180 <: AbstractOpticalElement end

"""
    Rotate180()

Optical element that computes the conjugate
field: E(x,y) → E*(x,y).

Physically, this is implemented by a phase conjugate mirror.
"""
struct Conjugate <: AbstractOpticalElement end

"""
    ConjugateInverter()

Optical element that computes the complex conjugate of the spatially
inverted field: E(r) → E*(-r).

In the coherence measurement setup, this is physically implemented by
a combination of a phase conjugate mirror and a beam inverter (Dove prism
or similar), or effectively by the geometry of the nonlinear crystal interaction.

Requires the field grid to be symmetric around zero.
"""
struct ConjugateInverter <: AbstractOpticalElement end

"""
    SLM(φ, θ_xy, θ_z)

Arbitrary phase mask defined by a matrix of phase values [rad].

Give the angles in degrees for the grating direction of the beam.
That is, how much does the first order of the grating differ from the zeroth order.
If θ_xy = θ_z = 0, then the zeroth order and first order point in the same direction.
"""
struct SLM <: AbstractOpticalElement
  φ :: Matrix{Float64}
  θ_xy :: Float64
  θ_z :: Float64
end

"""
    PhaseMask(φ)

Arbitrary phase mask defined by a matrix of phase values [rad].
"""
struct PhaseMask <: AbstractOpticalElement
  φ :: Matrix{Float64}
end
