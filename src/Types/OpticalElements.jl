"""
    FreeSpace(z)

Free space propagation by distance `z` [m].
"""
struct FreeSpace <: AbstractOpticalElement
    z :: Float64

    function FreeSpace(z::Float64)
        @assert z >= 0 "Propagation distance must be nonnegative"
        new(z)
    end
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

# TODO: PENDING... This can make trouble if we leave it as a matrix without units and 
# no way to match different grids.
"""
    PhaseMask(φ)

Arbitrary phase mask defined by a matrix of phase values [rad].
"""
struct PhaseMask <: AbstractOpticalElement
    φ :: Matrix{Float64}
end
