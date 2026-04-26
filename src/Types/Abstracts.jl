abstract type AbstractField end
abstract type AbstractBeam end
abstract type AbstractGrid end

grid(f::AbstractField) = error("$(typeof(f)) must implement grid()")
wavelength(f::AbstractField) = error("$(typeof(f)) must implement wavelength()")
intensity(f::AbstractField) = error("$(typeof(f)) must implement intensity()")

wavelength(b::AbstractBeam) = error("$(typeof(b)) must implement wavelength()")
evaluate(g::AbstractGrid, b::AbstractBeam) = error("$(typeof(b)) must implement sample()")

abstract type AbstractOpticalElement end
abstract type AbstractEnsembleSource <: AbstractBeam end

##########
# Checks
##########

function same_grid(f1::AbstractField, f2::AbstractField)
  return grid(f1) === grid(f2)
end

function same_wavelength(f1::AbstractField, f2::AbstractField)
  return wavelength(f1) ≈ wavelength(f2)
end

function check_compatible(f1::AbstractField, f2::AbstractField)
  same_grid(f1, f2) || error("Fields must share the same grid")
  same_wavelength(f1, f2) || error("Fields must have the same wavelength")
end


