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
include("Types/Mixers.jl")

include("Utils.jl")

# Implementations

include("Optics/Elements.jl")
include("Optics/Mixers.jl")
#
include("Sources/Ensembles.jl")
include("Sources/ShellModels.jl")
#
include("Propagation/Fresnel.jl") # Fresnel diffraction
#include("Propagation/Fraunhofer.jl") # Far-field / Fraunhofer diffraction
include("Propagation/AngularSpectrum.jl") # Angular spectrum method
include("Propagation/Bluestein.jl")
#include("Propagation/BeamPropagation.jl") # BPM, split-step methods
#
include("Beams/StructuredLight.jl") # LG, HG, Bessel, Airy beams
#
#include("Quantum/States.jl") 
#include("Quantum/Operators.jl") 
#include("Quantum/Cavity.jl") 
#
#include("Solvers/FDTD.jl")
#include("Solvers/Eigenmode.jl")
#include("Solvers/Optimization.jl")
#
include("Setups/Setups.jl")
#
include("Visualization/Recipes.jl") 
include("Visualization/Plots.jl") 
#include("Visualization/Animations.jl") 

# TODO: Add file for Bluestein propagation method and chirp z transform
# based on Hu et al. 2020 - Efficient Full-path optical calculation of scalar
# and vector diffraction using the Bluestein method.

export AbstractGrid, AbstractField, AbstractBeam
export AbstractOpticalElement, AbstractEnsembleSource, AbstractMixer

export TransverseGrid

export ScalarField

export SphericalBeam, ParaboloidalBeam, GaussianBeam, LGBeam, HGBeam

export grid, wavelength, intensity # Beams
export evaluate                                     # Grid, Beam -> ScalarField

export apply, propagate_angular, propagate_fresnel, propagate_bluestein

export EnsembleIntensity
export sample, ensemble_intensity  # Ensembles

export mix

export FreeSpace, ThinLens, FourierLens, CircularAperture
export SpiralPhaseElement, SLM, PhaseMask
export Rotate180, Conjugate, ConjugateInverter

export GaussianSchellModel, DiskSchellModel

export SFGCrystal

export OpticalSystem

export mixing_experiment

export get_wavevectors

end
