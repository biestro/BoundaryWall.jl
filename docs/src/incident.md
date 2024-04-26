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

Here, we deal with only real valued angle beam shaping (contrary to )