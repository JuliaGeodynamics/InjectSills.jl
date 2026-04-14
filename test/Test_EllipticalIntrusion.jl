using Test
using GeoParams, InjectSills

CharDim = GEO_units(length=1000m, temperature=1000C, stress=10Pa, viscosity=1e20Pas)

e2 = EllipticalIntrusion(Center=Point2(0.0, -5000.0)*m, Angle=Vec1(0.0)*NoUnits, W=2000.0m, H=200.0m)
@test isdimensional(e2) == true
@test isdimensional(nondimensionalize(e2, CharDim)) == false

e2u = EllipticalIntrusion(e2, W=2500.0m, H=250.0m)
@test UnitValue(e2u.W) ≈ 2500.0m
@test UnitValue(e2u.H) ≈ 250.0m

@test InjectSills.area(e2) ≈ π * UnitValue(e2.W) / 2 * UnitValue(e2.H) / 2
@test InjectSills.volume(e2) ≈ (4 / 3) * π * (UnitValue(e2.W) / 2)^2 * (UnitValue(e2.H) / 2)

d_center = hostrock_displacement(e2, Point2(0.0, -5000.0))
@test d_center[1] ≈ 0.0
@test d_center[2] ≈ 0.0

d_off = hostrock_displacement(e2, Point2(200.0, -4980.0))
@test !isnan(d_off[1])
@test !isnan(d_off[2])

@test inside(Point2(0.0, -5000.0), e2) == true
@test inside(Point2(1000.0, -5000.0), e2) == true
@test inside(Point2(1001.0, -5000.0), e2) == false
@test inside(Point2(0.0, -4900.0), e2) == true
@test inside(Point2(0.0, -4899.0), e2) == false

# Rotation should produce non-axis-aligned displacement components
e2r = EllipticalIntrusion(Center=Point2(0.0, -5000.0)*m, Angle=Vec1(45.0)*NoUnits, W=2000.0m, H=200.0m)
d_rot = hostrock_displacement(e2r, Point2(200.0, -5010.0))
@test abs(d_rot[1]) > 0.0
@test abs(d_rot[2]) > 0.0

e3 = EllipticalIntrusion(Center=Point3(0.0, 0.0, -5000.0)*m, Angle=Vec2(0.0, 0.0)*NoUnits, W=2000.0m, H=200.0m)
@test inside(Point3(0.0, 0.0, -5000.0), e3) == true
@test inside(Point3(1000.0, 0.0, -5000.0), e3) == true
@test inside(Point3(1001.0, 0.0, -5000.0), e3) == false

d3 = hostrock_displacement(e3, Point3(200.0, 100.0, -4980.0))
@test !any(isnan, Tuple(d3))
