"""
    mixing_experiment(grid, source, system1, system2, mixer, z;
                      n_realizations=500) -> Matrix{ComplexF64}

Monte Carlo simulation of a two-arm nonlinear mixing experiment using
a single spatially partially coherent source.

The source is sampled once per realization and the resulting field is
shared between both arms — physically corresponding to a beam splitter.
Each arm then propagates through its respective optical system before
being mixed in the nonlinear crystal.

For example, with `system2` containing a `ConjugateInverter` and an
`SFGCrystal` as mixer, the output is:

```math
\\langle E(\\mathbf{r}) E^*(-\\mathbf{r}) \\rangle = W(\\mathbf{r}, -\\mathbf{r})
```

# Arguments
- `grid`          : sampling grid at the source plane
- `source`        : the partially coherent source
- `system1`       : optical system for arm 1
- `system2`       : optical system for arm 2
- `system3`       : optical system for arm 3
- `mixer`         : nonlinear mixing element (e.g. `SFGCrystal`)
- `z`             : Distance to first element
- `n_realizations`: number of Monte Carlo realizations
"""
function mixing_experiment(
  grid           :: TransverseGrid{<:Real},
  source         :: AbstractEnsembleSource,
  system1        :: OpticalSystem,
  system2        :: OpticalSystem,
  system3        :: OpticalSystem,
  mixer          :: AbstractMixer,
  z              :: Float64;
  n_realizations :: Int,
) :: Tuple{ScalarField, ScalarField, ScalarField, Vector{Tuple{Float64, Float64}}}
  @assert n_realizations > 0

  grid_inv = TransverseGrid(-grid.x, -grid.y)

  W1 = zeros(ComplexF64, grid.Ny, grid.Nx)
  W2 = zeros(ComplexF64, grid.Ny, grid.Nx)
  W3 = zeros(ComplexF64, grid.Ny, grid.Nx)

  λ_in = source.λ
  λ_out = (λ_in^-1 + λ_in^-1)^-1

  points = Vector{Tuple{Float64,Float64}}(undef, n_realizations)

  for k in 1:n_realizations
    x0, y0 = _sample_disk(source.R)

    λ = source.λ
    beam1 = SphericalBeam(λ, z, center=(x0, y0))
    field1 = evaluate(grid, beam1)

    beam2 = SphericalBeam(λ, z, center=(-x0, -y0))
    field2 = evaluate(grid, beam2)

    ## Sample once — beam splitter gives same field to both arms
    #field, point = sample(source, grid, z)
    points[k] = (x0, y0)

    field1 = system1(field1)
    field2 = system2(field2)
    mixed = mix(mixer, field1, field2)
    field3 = system3(mixed)

    W1 += field1.U
    W2 += field2.U
    W3 += field3.U
  end

  field1 = ScalarField(W1 ./ n_realizations, grid, λ_in)
  field2 = ScalarField(W2 ./ n_realizations, grid, λ_in)
  field3 = ScalarField(W3 ./ n_realizations, grid, λ_out)

  return field1, field2, field3, points
end


"""
    mixing_experiment(grid, beam, system1, system2, mixer, z;
                      n_realizations=500) -> Matrix{ComplexF64}

Monte Carlo simulation of a two-arm nonlinear mixing experiment using
a single spatially partially coherent source.

The source is sampled once per realization and the resulting field is
shared between both arms — physically corresponding to a beam splitter.
Each arm then propagates through its respective optical system before
being mixed in the nonlinear crystal.

For example, with `system2` containing a `ConjugateInverter` and an
`SFGCrystal` as mixer, the output is:

```math
\\langle E(\\mathbf{r}) E^*(-\\mathbf{r}) \\rangle = W(\\mathbf{r}, -\\mathbf{r})
```

# Arguments
- `grid`          : sampling grid at the source plane
- `beamtype`      : the beam type
- `radius`        : disk source radius
- `system1`       : optical system for arm 1
- `system2`       : optical system for arm 2
- `system3`       : optical system for arm 3
- `mixer`         : nonlinear mixing element (e.g. `SFGCrystal`)
- `z`             : Distance to first element
- `n_realizations`: number of Monte Carlo realizations
"""
function mixing_experiment(
  grid           :: TransverseGrid{<:Real},
  lg1            :: LGBeam,
  lg2            :: LGBeam,
  radius         :: Float64,
  system1        :: OpticalSystem,
  system2        :: OpticalSystem,
  system3        :: OpticalSystem,
  mixer          :: AbstractMixer;
  n_realizations :: Int,
) :: Tuple{ScalarField, ScalarField, ScalarField, Vector{Tuple{Float64, Float64}}}
  @assert n_realizations > 0

  W1 = zeros(ComplexF64, grid.Ny, grid.Nx)
  W2 = zeros(ComplexF64, grid.Ny, grid.Nx)
  W3 = zeros(ComplexF64, grid.Ny, grid.Nx)

  λ1 = lg1.λ
  λ2 = lg2.λ
  λ_out = (λ1^-1 + λ2^-1)^-1

  points = Vector{Tuple{Float64,Float64}}(undef, n_realizations)

  for k in 1:n_realizations
    x0, y0 = _sample_disk(radius)
    points[k] = (x0, y0)

    beam1 = LGBeam(lg1.w0, lg1.λ, lg1.p, lg1.l, center=(x0, y0))
    field1 = evaluate(grid, beam1)

    beam2 = LGBeam(lg2.w0, lg2.λ, lg2.p, lg2.l, center=(x0, y0)) # Why don't I need -x0, -y0???
    field2 = evaluate(grid, beam2)

    field1 = system1(field1)
    field2 = system2(field2)
    mixed = mix(mixer, field1, field2)
    field3 = system3(mixed)

    W1 += field1.U
    W2 += field2.U
    W3 += field3.U
  end

  field1 = ScalarField(W1 ./ n_realizations, grid, λ1)
  field2 = ScalarField(W2 ./ n_realizations, grid, λ2)
  field3 = ScalarField(W3 ./ n_realizations, grid, λ_out)

  return field1, field2, field3, points
end
