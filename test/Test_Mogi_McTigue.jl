using Test
using GeoParams, InjectSills

CharDim = GEO_units(length=1000m, temperature=1000C, stress=10Pa, viscosity=1e20Pas)

# ---- MogiSphere ----------------------------------------------------------

mogi2D = MogiSphere(Center=Point2(0.0,-5000.0)*m, r=1500.0m, ΔP=10e6Pa, G=10e9Pa, ν=0.25*NoUnits)
@test isdimensional(mogi2D) == true

mogi2D_nd = nondimensionalize(mogi2D, CharDim)
@test isdimensional(mogi2D_nd) == false

# displacement at x=0 (directly above centre)
d = hostrock_displacement(mogi2D, Point2(0.0, 0.0))
@test d[1] ≈ 0.0            # no horizontal displacement on axis
@test d[2] ≈ 0.10125        # vertical (Uz), matches original Mogi_Uz[1]

# displacement at x=3500
d = hostrock_displacement(mogi2D, Point2(3500.0, 0.0))
@test d[2] ≈ 0.05566928318963278   # Uz, matches original Mogi_Uz[36]
@test d[1] ≈ 0.038968498232742954  # Ur, matches original Mogi_Ur[36]

# 3D sphere: axisymmetry means Ux==Uy on the x=y diagonal
mogi3D = MogiSphere(Center=Point3(0.0,0.0,-5000.0)*m, r=1500.0m, ΔP=10e6Pa, G=10e9Pa, ν=0.25*NoUnits)
d3 = hostrock_displacement(mogi3D, Point3(0.0, 0.0, 0.0))
@test d3[1] ≈ 0.0
@test d3[2] ≈ 0.0
@test d3[3] ≈ 0.10125

# array version
nx = 101
x = range(0.0, 10e3, length=nx)
Ux = zeros(nx); Uz = zeros(nx)
for i in eachindex(x)
    disp = hostrock_displacement(mogi2D, Point2(x[i], 0.0))
    Ux[i], Uz[i] = disp[1], disp[2]
end
@test Uz[1]  ≈ 0.10125
@test Ux[1]  ≈ 0.0

# inside
@test inside(Point2(0.0,   -5000.0), mogi2D) == true   # at centre
@test inside(Point2(1500.0,-5000.0), mogi2D) == true   # on boundary
@test inside(Point2(1501.0,-5000.0), mogi2D) == false  # just outside
@test inside(Point2(0.0,    0.0   ), mogi2D) == false  # at surface (far away)

@test inside(Point3(0.0, 0.0,   -5000.0), mogi3D) == true
@test inside(Point3(1500.0, 0.0,-5000.0), mogi3D) == true
@test inside(Point3(1501.0, 0.0,-5000.0), mogi3D) == false


# ---- McTigueSphere -------------------------------------------------------

mctigue2D = McTigueSphere(Center=Point2(0.0,-5000.0)*m, r=1500.0m, ΔP=10e6Pa, G=10e9Pa, ν=0.25*NoUnits)
@test isdimensional(mctigue2D) == true

mctigue2D_nd = nondimensionalize(mctigue2D, CharDim)
@test isdimensional(mctigue2D_nd) == false

# displacement at x=0
d = hostrock_displacement(mctigue2D, Point2(0.0, 0.0))
@test d[1] ≈ 0.0
@test d[2] ≈ 0.10407289402173914   # matches original McTigue_Uz[1]

# displacement at x=3500
d = hostrock_displacement(mctigue2D, Point2(3500.0, 0.0))
@test d[2] ≈ 0.05665722209549374   # McTigue_Uz[36]
@test d[1] ≈ 0.039660055466845624  # McTigue_Ur[36]

# inside (same geometry as Mogi)
@test inside(Point2(0.0,   -5000.0), mctigue2D) == true
@test inside(Point2(1500.0,-5000.0), mctigue2D) == true
@test inside(Point2(1501.0,-5000.0), mctigue2D) == false
