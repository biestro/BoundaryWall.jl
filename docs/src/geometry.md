# Geometry

## Defining curves

The effective use of the *Boundary Wall Method* (BWM) requires equal-length segments. 
and the need to interpolate uneven curves. While one might have an equispaced $x-y$ domain, the arc length $\mathrm{d}s$ need not be homogeneous.

The general steps for using this toolbox with a custom array of points are the following:
1. Define $(x_i,y_i)$
    - (Re)discretize using 

For closed curves, we implement the `discretizeCurve` function in order to interpolate the whole arc length into $N$ equal segments. For open curves, one must use `discretizeResonatorCurve` (it can be used for closed curves as well, but be mindful of using $N$ vs. $N+1$ points). Interpolating with `Interpolations.jl` will work as well.

Once one has a relatively accurate arc length discretization, one must calculate midpoints using the `calcMidpoints` function. Note that the total number of midpoints must strictly match the dimensions of the $\mathbf{T}$ matrix.

Finally, one might want to obtain the distance matrix representing $[\mathbf{r}]_{ij}=\mathbf{r}_{i}-\mathbf{r}_{j}$. The following example 

```@example midpoints
using BoundaryWall, CairoMakie # hide
θ     = LinRange(-pi, pi/2, 300)
x,y   = cos.(θ),sin.(θ)

x,y   = divideResonatorCurve(x,y, 300)
xm,ym = calcMidpoints(x, y)
ds    = calcArcLength(x,y)
rij   = calcDistances(xm, ym)
lines(ds)
```



A good way to check if one has done the discretization correctly is to check for the symmetry of the scattering $\mathbf{T}$ matrix, or to graph `diff(ds)` and check for homogeneous behaviour.

## Builtin methods

The package provides some useful functions such as `createEllipse` and `createCircle`. 