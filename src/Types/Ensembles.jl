"""
    EnsembleIntensity(I, I_std, grid, points)
"""
struct EnsembleIntensity
  I :: Matrix{Float64}
  I_std :: Matrix{Float64}
  grid :: TransverseGrid
  points :: Vector{Tuple{Float64,Float64}}
end
