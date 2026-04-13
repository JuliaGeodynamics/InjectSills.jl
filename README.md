[![CI](https://github.com/JuliaGeodynamics/InjectSills.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaGeodynamics/InjectSills.jl/actions/workflows/CI.yml)

# InjectSills.jl

This collects analytical solutions that give the displacement around sills, dikes, and magmatic intrusions embedded in an elastic half-space.
It can be combined with [GeophysicalModelGenerator](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl) to plot surface displacement caused by intrusion events, and with [JustPIC](https://github.com/JuliaGeodynamics/JustPIC.jl) to advect particles during sill injection.

## Available source types

All source types are subtypes of `AbstractSill{N,T}` and share a common interface: `hostrock_displacement`, `inside`, and `isdimensional`. Parameters are stored as [`GeoParams`](https://github.com/JuliaGeodynamics/GeoParams.jl) units and can be non-dimensionalized.

#### `PennyShapedSill` — penny-shaped tensile crack (2D / 3D)

Displacement field around a pressurized elliptical sill in a homogeneous elastic half-space.

```julia
sill = PennyShapedSill(
    Center = Point3(0.0, 0.0, -5000.0)*m,
    W      = 2000.0m,     # radius
    H      = 100.0m,      # maximum opening thickness
    ΔP     = 1e6Pa,
    E      = 1.5e10Pa,
    ν      = 0.3*NoUnits,
    Angle  = Vec2(10.0, 0.0),   # dip, strike [degrees]
)
d = hostrock_displacement(sill, Point3(x, y, z))
```

Reference: Sun, R.J. (1969): Theoretical size of hydraulically induced horizontal fractures and corresponding surface uplift in an idealized medium. *J. Geophys. Res.* 74, 5995–6011. https://doi.org/10.1029/JB074i025p05995

---

#### `MogiSphere` — Mogi point-source sphere (2D / 3D)

First-order solution for a pressurized spherical cavity; displacement is proportional to `(p − Center) / |p − Center|³`.

```julia
src = MogiSphere(
    Center = Point3(0.0, 0.0, -5000.0)*m,
    r      = 1500.0m,
    ΔP     = 10e6Pa,
    G      = 10e9Pa,
    ν      = 0.25*NoUnits,
)
d = hostrock_displacement(src, Point3(x, y, z))
```

Reference: Mogi, K. (1958): Relations between the eruptions of various volcanoes and the deformations of the ground surfaces around them. *Bull. Earthq. Res. Inst. Univ. Tokyo* 36, 99–134.

---

#### `McTigueSphere` — McTigue finite-radius sphere (2D / 3D)

Same geometry as `MogiSphere` with a higher-order `(r/d)³` correction for finite sphere size.

```julia
src = McTigueSphere(
    Center = Point3(0.0, 0.0, -5000.0)*m,
    r      = 1500.0m,
    ΔP     = 10e6Pa,
    G      = 10e9Pa,
    ν      = 0.25*NoUnits,
)
d = hostrock_displacement(src, Point3(x, y, z))
```

Reference: McTigue, D.F. (1987): Elastic stress and deformation near a finite spherical magma body: resolution of the point source paradox. *J. Geophys. Res.* 92 (B12), 12931–12940. https://doi.org/10.1029/JB092iB12p12931

---

#### `FiniteEllipsoidalCavity` — finite ellipsoidal cavity (3D)

Surface deformation from a uniformly-pressurized triaxial ellipsoidal cavity. Uses an adaptive grid of point Compound Dislocation Models (pCDMs) internally. The array form is efficient (grid built once per struct).

```julia
fec = FiniteEllipsoidalCavity(
    Center = Point3(0.0, 0.0, -10000.0)*m,
    ax     = 500.0m,  ay = 500.0m,  az = 2000.0m,
    Angle  = Vec{3}(0.0, 25.0, -10.0)*NoUnits,  # omegaX, omegaY, omegaZ [degrees]
    ΔP     = 10e6Pa,
    mu     = 10e9Pa,
    lambda = 10e9Pa,
)
# Array of surface observation points (efficient):
ue, un, uv, dV, DV, Ns = hostrock_displacement(fec, X, Y)

# Single point:
d = hostrock_displacement(fec, Point3(x, y, 0.0))
```

References:
- Nikkhoo, M., Rivalta, E. (2023): Surface deformations and gravity changes caused by pressurized finite ellipsoidal cavities. *Geophys. J. Int.* https://doi.org/10.1093/gji/ggac351
- Nikkhoo, M., Rivalta, E. (2022): Analytical solutions for gravity changes caused by triaxial volumetric sources. *Geophys. Res. Lett.* https://doi.org/10.1029/2021GL095442
- Nikkhoo, M., Walter, T.R., Lundgren, P.R., Prats-Iraola, P. (2017): Compound dislocation models (CDMs) for volcano deformation analyses. *Geophys. J. Int.* https://doi.org/10.1093/gji/ggw427

---

## Common interface

```julia
# Check if a point is inside the source volume
inside(Point3(x, y, z), src)

# Non-dimensionalize all parameters
src_nd = nondimensionalize(src, CharDim)
```

## `hostrock_displacement` — computing displacement fields

For all source types, `hostrock_displacement` can be called on a single point or on arrays.

**2D — single point**
```julia
Ux, Uz = hostrock_displacement(src, Point2(x, z))
```

**2D — array**
```julia
X = range(-10e3, 10e3, length=128)
Z = range(-20e3,    0, length=128)
Ux, Uz = hostrock_displacement(src, X, Z)   # returns two arrays
```

**3D — single point**
```julia
Ux, Uy, Uz = hostrock_displacement(src, Point3(x, y, z))
```

**3D — array**
```julia
X = range(-10e3, 10e3, length=64)
Y = range(-10e3, 10e3, length=64)
Z = range(-20e3,    0, length=64)
Ux, Uy, Uz = hostrock_displacement(src, X, Y, Z)   # returns three arrays
```

For `FiniteEllipsoidalCavity` the array form additionally returns the volume change, potency, and number of point CDMs:
```julia
ue, un, uv, dV, DV, Ns = hostrock_displacement(fec, X, Y)
```

## `surface_displacement` — displacement on a GeophysicalModelGenerator surface

When [GeophysicalModelGenerator](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl) is loaded, `surface_displacement` computes the displacement at every point of a `CartData` surface. `CartData` uses km; the conversion to metres is handled internally.

```julia
using InjectSills, GeophysicalModelGenerator

# Build a CartData surface (flat, at sea level)
nx, ny, nz = 201, 201, 1
x2D = [xi for xi in range(-100.0, 100.0, length=nx), _ in 1:ny,  _ in 1:nz]  # km
y2D = [yi for _ in 1:nx, yi in range(-100.0, 100.0, length=ny),  _ in 1:nz]
z2D = zeros(nx, ny, nz)   # km; replace with topography if available
surf = CartData(x2D, y2D, z2D, (Elevation = z2D,))

# Define a source
src = MogiSphere(
    Center = Point3(0.0, 0.0, -5000.0)*m,
    r = 1500.0m, ΔP = 10e6Pa, G = 10e9Pa, ν = 0.25*NoUnits,
)

# Option 1 — return displacement arrays (metres), useful for inversion
Ux, Uy, Uz = surface_displacement(src, surf)

# Option 2 — add displacement as a vector field to the CartData
surf_out = surface_displacement(src, surf; add_fields = true)
# surf_out.fields.Displacement_m = (Ux, Uy, Uz) tuple in metres
# components: [1] = Ux, [2] = Uy, [3] = Uz
```

The `add_fields = true` form preserves all existing fields in `surf` and appends `:Displacement_m` as a `(Ux, Uy, Uz)` tuple — consistent with the GMG convention for vector fields. The plain array form is preferable for inversion workflows where you want to avoid constructing new `CartData` objects at each iteration.
