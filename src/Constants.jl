module Constants

# Mathematical
const π = Base.π
const τ = 2π

# Fundamental physics
const c = 2.99792458e8 # speed of light [m/s] (exact)
const e = 1.602176634e-19 # elementary charge [C] (exact)
const kB = 1.380649e-23 # Boltzmann constant [J/K] (exact)
const ε₀ = 8.8541878128e-12 # vacuum permittivity [F/m]
const μ₀ = 1.25663706212e-6 # vacuum permeability [H/m]

# Optics
const λ_vis_min = 380e-9 # Visible spectrum min [m]
const λ_vis_max = 780e-9 # Visible spectrum max [m]

# Common laser wavelengths [m]
const λ_HeNe = 632.8e-9 # HeNe laser
const λ_Nd_YAG = 1064e-9 # Nd:YAG fundamental
const λ_Nd_YAG2 = 532e-9 # Nd:YAG second harmonic
const λ_Ti_Sa_c = 800e-9 # Ti:Sapphire
const λ_blue = 405e-9 # Blue-ray / GaN diode
const λ_CO2 = 10.6e-6 # CO₂ laser

# Common refractive indexes at 589 nm (sodium D line)
const n_vacuum = 1.0 # vacuum
const n_air = 1.0003 # refractive index of air at STP
const n_silica = 1.4585 # fused silica
const n_BK7 = 1.5168 # Schott BK7
const n_sapphire = 1.7682 # Al₂O₃

end
