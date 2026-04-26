# Examples

## Gaussian Beams
```@example GaussianBeam
using Optikos
using CairoMakie

grid = TransverseGrid(range(-1e-3, 1e-3, 256), range(-1e-3, 1e-3, 256))
beam = GaussianBeam(200e-6, 632.8e-9)
field = evaluate(grid, beam)

fig, ax, hm = fieldplot(grid, beam, units=:mm)
ax.xlabel = "x [mm]"
ax.ylabel = "y [mm]"
ax.title  = "Gaussian beam intensity"
Colorbar(fig[1,2], hm)
fig
```

## Laguerre-Gaussian Beams
```@example LGBeam
using Optikos
using CairoMakie

grid = TransverseGrid(range(-1e-3, 1e-3, 256), range(-1e-3, 1e-3, 256))
beam = LGBeam(200e-6, 632.8e-9, 0, 1)
field = evaluate(grid, beam)

fig, ax, hm = fieldplot(field, units=:mm)
ax.xlabel = "x [mm]"
ax.ylabel = "y [mm]"
ax.title  = "LG beam (p=0, l=1) intensity"
Colorbar(fig[1,2], hm)
fig
```
