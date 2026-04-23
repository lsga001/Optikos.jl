module Optikos

# Basics
# TODO: Create a proper project structure

include("Constants.jl")
include("Types.jl")
include("Utils.jl")

#include("Optics/General.jl")
#include("Optics/Ray.jl")
#include("Optics/Wave.jl")
#
#include("Propagation/Fresnel.jl") # Fresnel diffraction
#include("Propagation/Fraunhofer.jl") # Far-field / Fraunhofer diffraction
#include("Propagation/AngularSpectrum.jl") # Angular spectrum method
#include("Propagation/BeamPropagation.jl") # BPM, split-step methods
#
#include("Beams/Gaussian.jl") # Gaussian beams, beam waist, Rayleigh range
#include("Beams/Structured.jl") # LG, HG, Bessel, Airy beams
#include("Beams/Polarization.jl") # Stokes params, Jones calculus, Mueller matrices
#
#include("Quantum/States.jl") 
#include("Quantum/Operators.jl") 
#include("Quantum/Cavity.jl") 
#
#include("Solvers/FDTD.jl")
#include("Solvers/Eigenmode.jl")
#include("Solvers/Optimization.jl")
#
#include("Visualization/Plots.jl") 
#include("Visualization/Animations.jl") 

function operateOnXY(x=1, y=2)
  println("This function is modified")
  return x*y
end

end
