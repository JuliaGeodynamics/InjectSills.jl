using Test
using GeoParams, InjectSills

CharDim = GEO_units(length=1000m, temperature=1000C, stress=10Pa, viscosity=1e20Pas)

s2 = SquareDikeTopAccretion(Center=Point2(0.0, -5000.0)*m, Angle=Vec1(0.0)*NoUnits, W=2000.0m, H=100.0m)
@test isdimensional(s2) == true

s2_nd = nondimensionalize(s2, CharDim)
@test isdimensional(s2_nd) == false

s2_u = SquareDikeTopAccretion(s2, W=2500.0m, H=120.0m)
@test UnitValue(s2_u.W) ≈ 2500.0m
@test UnitValue(s2_u.H) ≈ 120.0m

s2_u2 = SquareDikeTopAccretion(s2, H=80.0)
@test UnitValue(s2_u2.H) ≈ 80.0m

@test InjectSills.area(s2) ≈ UnitValue(s2.W) * UnitValue(s2.H)
@test InjectSills.volume(s2) ≈ UnitValue(s2.W)^2 * UnitValue(s2.H)

# only lower half displaces downwards
d_below = hostrock_displacement(s2, Point2(0.0, -5000.0 - 10.0))
@test d_below[1] ≈ 0.0
@test d_below[2] ≈ -100.0

d_above = hostrock_displacement(s2, Point2(0.0, -5000.0 + 10.0))
@test d_above[1] ≈ 0.0
@test d_above[2] ≈ 0.0

d_out = hostrock_displacement(s2, Point2(2000.0, -5000.0))
@test d_out[1] ≈ 0.0
@test d_out[2] ≈ 0.0

s2_rot = SquareDikeTopAccretion(Center=Point2(0.0, -5000.0)*m, Angle=Vec1(45.0)*NoUnits, W=2000.0m, H=100.0m)
d_rot = hostrock_displacement(s2_rot, Point2(0.0, -5010.0))
@test abs(d_rot[1]) > 1.0
@test abs(d_rot[2]) > 1.0

@test inside(Point2(0.0, -5000.0), s2) == true
@test inside(Point2(999.0, -5000.0), s2) == true
@test inside(Point2(1001.0, -5000.0), s2) == false
@test inside(Point2(0.0, -4950.0), s2) == true
@test inside(Point2(0.0, -4949.0), s2) == false

s3 = SquareDikeTopAccretion(Center=Point3(0.0, 0.0, -5000.0)*m, Angle=Vec2(0.0, 0.0)*NoUnits, W=2000.0m, H=100.0m)
@test inside(Point3(0.0, 0.0, -5000.0), s3) == true
@test inside(Point3(999.0, 999.0, -5000.0), s3) == true
@test inside(Point3(1001.0, 0.0, -5000.0), s3) == false
@test inside(Point3(0.0, 0.0, -4950.0), s3) == true
@test inside(Point3(0.0, 0.0, -4949.0), s3) == false

d3_below = hostrock_displacement(s3, Point3(0.0, 0.0, -5010.0))
@test d3_below[1] ≈ 0.0
@test d3_below[2] ≈ 0.0
@test d3_below[3] ≈ -100.0

d3_above = hostrock_displacement(s3, Point3(0.0, 0.0, -4990.0))
@test d3_above[1] ≈ 0.0
@test d3_above[2] ≈ 0.0
@test d3_above[3] ≈ 0.0
