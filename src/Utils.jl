function get_wavevectors(k, θ_xy, θ_z)
  kx = k*sind(θ_z)*cosd(θ_xy)
  ky = k*sind(θ_z)*sind(θ_xy)
  kz = k*cosd(θ_z)
  return kx, ky, kz
end
