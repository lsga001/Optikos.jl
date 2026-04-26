module Optikos

using FFTW

# TODO: Create a proper project structure

# Basics
include("Constants.jl")

include("Types/Abstracts.jl")
include("Types/Grids.jl")
include("Types/Fields.jl")
include("Types/Beams.jl")
include("Types/Sources.jl")
include("Types/OpticalElements.jl")
include("Types/OpticalSystems.jl")
include("Types/Ensembles.jl")

include("Utils.jl")

# Implementations

include("Optics/Elements.jl")
#include("Optics/General.jl")
#include("Optics/Ray.jl")
#include("Optics/Wave.jl")
#include("Optics/Fourier.jl")
#include("Optics/Nonlinear.jl")
#
include("Sources/Ensembles.jl")
include("Sources/ShellModels.jl")
#
include("Propagation/Fresnel.jl") # Fresnel diffraction
#include("Propagation/Fraunhofer.jl") # Far-field / Fraunhofer diffraction
include("Propagation/AngularSpectrum.jl") # Angular spectrum method
#include("Propagation/BeamPropagation.jl") # BPM, split-step methods
#
#include("Beams/Gaussian.jl") # Gaussian beams, beam waist, Rayleigh range
include("Beams/StructuredLight.jl") # LG, HG, Bessel, Airy beams
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
include("Visualization/Recipes.jl") 
#include("Visualization/Plots.jl") 
#include("Visualization/Animations.jl") 

# TODO: Add file for Bluestein propagation method and chirp z transform
# based on Hu et al. 2020 - Efficient Full-path optical calculation of scalar
# and vector diffraction using the Bluestein method.

export AbstractGrid, AbstractField, AbstractBeam, AbstractOpticalElement, AbstractEnsembleSource
export grid, wavelength, intensity        # Beams
export evaluate                           # Grid, Beam -> ScalarField

export apply, propagate_angular, propagate_fresnel

export EnsembleIntensity
export sample, ensemble_intensity  # Ensembles

export TransverseGrid
export ScalarField
export SphericalBeam, ParaboloidalBeam, GaussianBeam, LGBeam, HGBeam
export FreeSpace, ThinLens, CircularAperture, PhaseMask
export GaussianSchellModel, DiskSchellModel
export OpticalSystem

end
