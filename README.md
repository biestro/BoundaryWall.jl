# BoundaryWall.jl

Simple toolbox that implements an integral method for solving for an incident (TM) wave with well behaved scatterers. Documentation (some of it at least) available at [https://biestro.github.io/posts/boundary-wall/](https://biestro.github.io/posts/boundary-wall/).

## Overview

**BoundaryWall.jl** provides a collection of algorithms and tools for 
calculating the scattering of arbitrary boundaries, based upon the work of
[M. G. E. da Luz, A. S. Lupu-Sax, and E. J. Heller on Quantum Scattering](https://doi.org/10.1103/PhysRevE.56.2496), 
written in pure Julia.

The design of this project was motivated for calculating eigenstates 
and scattering wavefunctions in billiards (see [F.M. Zanetti et al.](https://doi.org/10.1016/j.aop.2008.01.008) for original algorithm, our work in 
[A. Ruiz-Biestro, JC Gutiérrez-Vega](https://doi.org/10.1103/PhysRevE.109.034203)).

This raises an important detail about possible applications. This is a **Boundary Integral Method**,
which limits the applicable domains to 2D curves consisting of: [simply 
connected curves](https://mathworld.wolfram.com/SimplyConnected.html), open curves, or a collection of the previous two. 
Boundaries with intersections and inner loops **will generally not work**.
While the method does not require a parametrization of the boundary, it
does require certain homogeneity in arc length. For more details about
the implementation and theoretical foundations, see the [Geometry](geometry.md)
section.

![Permissible shapes](docs/src/assets/drawing.svg)

If you have questions or would like to brainstorm ideas, don't hesitate 
to contact me. I'd be happy to talk about the project and improve the
toolbox in general.

![Eigenstate](docs/src/assets/parabolic.gif)

## Installation

```julia
]add https://github.com/ModifiedBear/BoundaryWall.jl.git
using BoundaryWall 
```

## Documentation

Tutorials and documentation can be found in [biestro.github.io/software](https://biestro.github.io/software)


# TODO
- Fix documentation
