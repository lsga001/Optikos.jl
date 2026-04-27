"""
    ensemble_average(source, grid; n_realizations=500,
                     system=identity, compute_std=false) -> EnsembleIntensity

Compute the incoherently averaged intensity of a partially coherent
source by Monte Carlo ensemble averaging.

Each realization is:
1. Drawn from `source` via `sample`
2. Passed through `system`
3. Its intensity |E|² is accumulated via Welford's algorithm

# Arguments
- `source`        : an `AbstractEnsembleSource`
- `grid`          : `TransverseGrid` for field sampling
- `z`             : Initial sampling distance
- `n_realizations`: number of Monte Carlo realizations
- `system`        : any callable `ScalarField → ScalarField`.
                    Defaults to `identity`.
- `compute_std`   : if `true`, also compute the standard deviation
                    of intensity across realizations

# References
Goodman, *Statistical Optics*, 2nd ed., §5.2
"""
function ensemble_intensity(
    source        :: AbstractEnsembleSource,
    grid          :: TransverseGrid{<:Real},
    z             :: Float64;
    n_realizations:: Int  = 500,
    system              = identity,
    compute_std   :: Bool = false,
) :: EnsembleIntensity

  @assert n_realizations > 1 "Need at least 2 realizations"

  mean_I = zeros(Float64, grid.Ny, grid.Nx)
  M2     = compute_std ? zeros(Float64, grid.Ny, grid.Nx) : nothing
  points = Vector{Tuple{Float64,Float64}}(undef, n_realizations)

  for k in 1:n_realizations
    field, point = sample(source, grid, z)
    points[k] = point
    field = system(field)
    I_k   = intensity(field)

    # Welford's online algorithm
    δ      = I_k .- mean_I
    mean_I .+= δ ./ k
    if compute_std
        M2 .+= δ .* (I_k .- mean_I)
    end
  end

  std = compute_std ?
    sqrt.(M2 ./ (n_realizations - 1)) :
    nothing

  return EnsembleIntensity(mean_I, std, grid, points)
end

