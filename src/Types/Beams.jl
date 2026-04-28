"""
    SphericalBeam(λ::Float64, z; n_index=1.0, center=(0.0, 0.0))

Analytical descriptor of a spherical beam.

# Fields
- `λ `: wavelength [m]
- `z`: axial position of the transverse beam profile assuming origin at z=0 [m]
- `n_index `: refractive index of the propagation medium

# Examples
```jldoctest
julia> beam = SphericalBeam(632.8e-9, 0.1);
```

# References
Saleh & Teich, *Fundamentals of Photonics*, 3rd ed., §2.2
"""
struct SphericalBeam <: AbstractBeam
  λ   ::Float64
  z  ::Float64
  n_index   ::Float64
  center :: Tuple{Float64, Float64}
end

function SphericalBeam(λ, z; n_index=1.0, center=(0.0, 0.0))
  @assert λ > 0 "Wavelength must be positive"
  @assert Base.abs(z) > 0 "Plane z-position must be nonzero"
  @assert n_index > 0 "Refractive index must be positive"
  return SphericalBeam(λ, z, n_index, center)
end

wavelength(beam::SphericalBeam) = beam.λ

"""
    ParaboloidalBeam(λ::Float64, z; n_index=1.0, center=(0.0, 0.0))

Analytical descriptor of a paraboloidal beam.

# Fields
- `λ `: wavelength [m]
- `z`: axial position of the transverse beam profile assuming origin at z=0 [m]
- `n_index `: refractive index of the propagation medium

# Examples
```jldoctest
julia> beam = ParaboloidalBeam(632.8e-9, 0.1);
```

# References
Saleh & Teich, *Fundamentals of Photonics*, 3rd ed., §2.2
"""
struct ParaboloidalBeam <: AbstractBeam
  λ         ::Float64
  z         ::Float64
  n_index   ::Float64
  center :: Tuple{Float64, Float64}
end

function ParaboloidalBeam(λ, z; n_index=1.0, center=(0.0, 0.0))
  @assert λ > 0 "Wavelength must be positive"
  @assert Base.abs(z) > 0 "Plane z-position must be nonzero"
  @assert n_index > 0 "Refractive index must be positive"
  return ParaboloidalBeam(λ, z, n_index, center)
end

wavelength(beam::ParaboloidalBeam) = beam.λ

"""
    GaussianBeam(w0::Float64, λ::Float64; z0=0.0, n_index=1.0, center=(0.0,0.0))

Analytical descriptor of a fundamental Gaussian beam (TEM₀₀).

# Fields
- `w0`: beam waist radius [m]
- `λ `: wavelength [m]
- `z0`: axial position of the transverse beam profile assuming origin at z=0 [m]
- `n_index `: refractive index of the propagation medium

# Examples
```jldoctest
julia> beam = GaussianBeam(200e-6, 632.8e-9);

julia> beam.w0
0.0002

julia> beam.n_index
1.0
```

# References
Saleh & Teich, *Fundamentals of Photonics*, 3rd ed., §3.1
"""
struct GaussianBeam <: AbstractBeam
  w0  ::Float64
  λ   ::Float64
  z0  ::Float64
  n_index   ::Float64
  center :: Tuple{Float64, Float64}

end

function GaussianBeam(w0, λ; z0=0.0, n_index=1.0, center=(0.0, 0.0))
  @assert w0 > 0 "Beam waist must be positive"
  @assert λ > 0 "Wavelength must be positive"
  @assert n_index > 0 "Refractive index must be positive"
  return GaussianBeam(w0, λ, z0, n_index, center)
end

wavelength(beam::GaussianBeam) = beam.λ


"""
    LGBeam(w0::Float64, λ::Float64, p::Int, l::Int; z0=0.0, n_index=1.0, center=(0.0, 0.0))

Analytical descriptor of a Laguerre-Gaussian beam LG(p, l).

Carries orbital angular momentum of `l·ħ` per photon. The azimuthal
index `l` determines the topological charge of the optical vortex,
and `p` is the number of radial nodes.

# Fields
- `w0`: beam waist radius [m]
- `λ `: wavelength [m]
- `p `: radial index (p ≥ 0)
- `l `: azimuthal index / topological charge (integer)
- `z0`: axial position of the transverse beam profile assuming origin at z=0 [m]
- `n_index `: refractive index of the propagation medium

# Examples
```jldoctest
julia> beam = LGBeam(100e-6, 532e-9, 0, 1);

julia> beam.l   # topological charge
1

julia> beam.p   # radial index
0
```

# References
Allen et al., *Phys. Rev. A* 45, 8185 (1992)
Saleh & Teich, *Fundamentals of Photonics*, 3rd ed., §3.3
"""
struct LGBeam <: AbstractBeam
  w0  ::Float64
  λ   ::Float64
  p   ::Int
  l   ::Int
  z0  ::Float64
  n_index   ::Float64
  center :: Tuple{Float64, Float64}

end

function LGBeam(w0, λ, p, l; z0=0.0, n_index=1.0, center=(0.0, 0.0))
  @assert w0 > 0 "Beam waist must be positive"
  @assert λ  > 0 "Wavelength must be positive"
  @assert p  >= 0 "Radial index p must be non-negative"
  @assert n_index  > 0 "Refractive index must be positive"
  return LGBeam(w0, λ, p, l, z0, n_index, center)
end

wavelength(beam::LGBeam) = beam.λ


"""
    HGBeam(w0::Float64, λ::Float64, l::Int, m::Int; z0=0.0, n_index=1.0, center=(0.0, 0.0))

Analytical descriptor of a Hermite-Gaussian beam HG(l, m).

Rectangular-symmetry transverse modes. The indices `l` and `m`
give the number of nodes along the x and y axes respectively.

# Fields
- `w0`   : beam waist radius [m]
- `λ`    : wavelength [m]
- `l`    : mode index along x (m ≥ 0)
- `m`    : mode index along y (n ≥ 0)
- `z0`: axial position of the transverse beam profile assuming origin at z=0 [m]
- `n_index` : refractive index of the propagation medium

# Examples
```jldoctest
julia> beam = HGBeam(150e-6, 1064e-9, 1, 0);

julia> beam.l
1

julia> beam.m
0
```

# References
Saleh & Teich, *Fundamentals of Photonics*, 3rd ed., §3.3
"""
struct HGBeam <: AbstractBeam
  w0  ::Float64
  λ   ::Float64
  l   ::Int
  m   ::Int
  z0  ::Float64
  n_index ::Float64
  center :: Tuple{Float64, Float64}
end

function HGBeam(w0, λ, l, m; z0=0.0, n_index=1.0, center=(0.0, 0.0))
  @assert w0 > 0  "Beam waist must be positive"
  @assert λ > 0  "Wavelength must be positive"
  @assert l >= 0 "Mode index m must be non-negative"
  @assert m >= 0 "Mode index n must be non-negative"
  @assert n_index > 0  "Refractive index must be positive"
  return HGBeam(w0, λ, l, m, z0, n_index, center)
end

wavelength(beam::HGBeam) = beam.λ
