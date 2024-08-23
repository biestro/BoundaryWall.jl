# Incident waves

Since the BWM solves the Lippmann-Schwinger Equation, which is essentially 
just an integral version of the Helmholtz equation, any solution to both
can be used as the incident field. Furthermore, any linear combination 
can also be used as the impinging wave.

```@docs
planeWave
# planeWave(k::SVector{2,Float64}, r::SVector{2,Float64})
```

## Beam shaping

Here, we deal with only real valued angle beam shaping, applied as a Gaussian beam.
Since a *plane wave* is a solution to the Helmholtz equation, therefore any linear
combination of these will also be a solution. Taken from
[Willis, et al. (2008)](https://opg.optica.org/oe/fulltext.cfm?uri=oe-16-3-1903&id=149977), one effectively *shapes* plane waves $\Phi(\mathbf{k}, \mathbf{r})$ with different angles and intensities into a Gaussian profile of width $\omega$,

```math
\Psi_z(\mathbf{r})=\frac{1}{2\pi}\int\mathrm{d}\mathbf{k} \; A(\mathbf{k})\,\Phi(\mathbf{k},\mathbf{r}),
```

where $A(\mathbf{k})$ is the spectra (a Gaussian). Note that, under this formalism, 
any other beam reconstructible from plane waves can be used as an incident wave.

```@docs
gaussianWave
```

!!! note
    
    Setting the precision for the quadrature is quite arbitrary, but necessary. 
    Since this is a Fourier decomposition, it will attempt to force periodicity
    at infinity, which implies a slow convergence rate.

```@example beam
using WGLMakie # hide
using BoundaryWall # hide
using StaticArrays: SVector # hide

x = LinRange(-5,5,200); # hide
y = LinRange(-2.5,2.5, 100); # hide
w = [gaussianWave(SVector(10.0, 0.0), SVector(_x,_y), 8.0; abstol=1e-7) for _x in x, _y in y]
fig = Figure() # hide
ax = Axis(fig[1,1]) # hide
heatmap!(ax, x, y, real(w), colormap=:balance) # hide
ax.aspect=DataAspect() # hide
fig # hide
```
