using Test
using GeoParams, InjectSills

CharDim = GEO_units(length=1000m, temperature=1000C, stress=10Pa, viscosity=1e20Pas)

# make mesh
x = 0:100:20e3
y = 0:200:20e3
X, Y = InjectSills.meshgrid(x, y)

# ---- construct and basic properties ------------------------------------

fec = FiniteEllipsoidalCavity(
    Center = Point3(0.0, 0.0, -10250.0)*m,
    ax     = 550.0m,
    ay     = 550.0m,
    az     = 3750.0m,
    Angle  = Vec{3}(0.0, 0.0, 0.0)*NoUnits,
    ΔP     = 18.8e6Pa,
    mu     = 10e9Pa,
    lambda = 10e9Pa,
    Nmax   = 5000,
    Cr     = 14,
)
@test isdimensional(fec) == true

fec_nd = nondimensionalize(fec, CharDim)
@test isdimensional(fec_nd) == false

# geometric consistency: 2D footprint area and 3D cavity volume
@test InjectSills.area(fec) ≈ π * UnitValue(fec.ax) * UnitValue(fec.ay)
@test InjectSills.volume(fec) ≈ (4 / 3) * π * UnitValue(fec.ax) * UnitValue(fec.ay) * UnitValue(fec.az)

# copy-constructor style updates
fec_u = FiniteEllipsoidalCavity(fec, ax=600.0m, ΔP=20e6Pa, Nmax=4000)
@test UnitValue(fec_u.ax) ≈ 600.0m
@test UnitValue(fec_u.ΔP) ≈ 20e6Pa
@test fec_u.Nmax == 4000

fec_u2 = FiniteEllipsoidalCavity(fec, ay=700.0)
@test UnitValue(fec_u2.ay) ≈ 700.0m
@test InjectSills.area(fec_u2) ≈ π * UnitValue(fec_u2.ax) * UnitValue(fec_u2.ay)
@test InjectSills.volume(fec_u2) ≈ (4 / 3) * π * UnitValue(fec_u2.ax) * UnitValue(fec_u2.ay) * UnitValue(fec_u2.az)

# ---- array displacement (matches original fECM results) ----------------

ue, un, uv, dV, DV, Ns = hostrock_displacement(fec, X, Y)

@test all(isapprox.(extrema(uv), (8.044701021001343e-04, 0.005989663901181), rtol=1e-4))
@test uv[1,1]      ≈ 0.004201484984743  rtol=1e-4
@test uv[9,100]    ≈ 0.004833273685574  rtol=1e-4
@test uv[87,46]    ≈ 0.002210474017935  rtol=1e-4
@test maximum(ue)  ≈ 0.005037030511448  rtol=1e-4
@test maximum(un)  ≈ 0.005037030511448  rtol=1e-4

# ---- single-point displacement -----------------------------------------

d = hostrock_displacement(fec, Point3(0.0, 0.0, 0.0))
@test d[3] ≈ uv[1,1]  rtol=1e-4   # vertical at (0,0) should match grid value

# ---- rotated cavity ----------------------------------------------------

fec2 = FiniteEllipsoidalCavity(
    Center = Point3(0.0, 0.0, -5000.0)*m,
    ax     = 850.0m,
    ay     = 550.0m,
    az     = 1500.0m,
    Angle  = Vec{3}(0.0, 25.0, -65.0)*NoUnits,
    ΔP     = 18.8e6Pa,
    mu     = 10e9Pa,
    lambda = 10e9Pa,
)
ue2, un2, uv2, _, _, _ = hostrock_displacement(fec2, X, Y)

@test all(isapprox.(extrema(uv2), (3.405925303200454e-04, 0.022797797325954), rtol=1e-4))
@test uv2[1,1]      ≈ 0.020284896020590  rtol=1e-4
@test uv2[9,100]    ≈ 0.004740361583832  rtol=1e-4
@test uv2[87,46]    ≈ 0.001025808009311  rtol=1e-4
@test maximum(ue2)  ≈ 0.017285000082854  rtol=1e-4
@test maximum(un2)  ≈ 0.011987093395042  rtol=1e-4

# ---- inside ------------------------------------------------------------

# at centre (no rotation)
fec3 = FiniteEllipsoidalCavity(
    Center = Point3(0.0, 0.0, -10000.0)*m,
    ax=500.0m, ay=500.0m, az=2000.0m,
    Angle=Vec{3}(0.0,0.0,0.0)*NoUnits,
    ΔP=1e6Pa, mu=10e9Pa, lambda=10e9Pa,
)
@test inside(Point3(0.0, 0.0, -10000.0), fec3) == true   # centre
@test inside(Point3(500.0, 0.0, -10000.0), fec3) == true  # on x boundary
@test inside(Point3(501.0, 0.0, -10000.0), fec3) == false # just outside x
@test inside(Point3(0.0, 0.0, -10000.0 + 2000.0), fec3) == true  # on z boundary
@test inside(Point3(0.0, 0.0, -10000.0 + 2001.0), fec3) == false # just outside z
@test inside(Point3(0.0, 0.0, 0.0), fec3) == false        # at surface (far above)
