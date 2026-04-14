using Test
using GeoParams, InjectSills

CharDim = GEO_units(length=1000m, temperature=1000C, stress=10Pa, viscosity=1e20Pas)

# ---- CylindricalDikeTopAccretion -----------------------------------------

c2 = CylindricalDikeTopAccretion(Center=Point2(0.0, -5000.0)*m, Angle=Vec1(0.0)*NoUnits, W=2000.0m, H=100.0m)
@test isdimensional(c2) == true
@test isdimensional(nondimensionalize(c2, CharDim)) == false

c2u = CylindricalDikeTopAccretion(c2, W=2500.0m, H=120.0m)
@test UnitValue(c2u.W) ≈ 2500.0m
@test UnitValue(c2u.H) ≈ 120.0m

@test InjectSills.area(c2) ≈ UnitValue(c2.W) * UnitValue(c2.H)
@test InjectSills.volume(c2) ≈ π * (UnitValue(c2.W) / 2)^2 * UnitValue(c2.H)

d_below = hostrock_displacement(c2, Point2(100.0, -5010.0))
@test d_below[1] ≈ 0.0
@test d_below[2] ≈ -100.0

d_above = hostrock_displacement(c2, Point2(100.0, -4990.0))
@test d_above[2] ≈ 0.0

# MTK behavior: condition is x <= W/2 (not abs(x)), so far negative x is also "inside"
d_left = hostrock_displacement(c2, Point2(-5000.0, -5010.0))
@test d_left[2] ≈ -100.0

@test inside(Point2(0.0, -5000.0), c2) == true
@test inside(Point2(999.0, -5000.0), c2) == true
@test inside(Point2(1001.0, -5000.0), c2) == false

c3 = CylindricalDikeTopAccretion(Center=Point3(0.0, 0.0, -5000.0)*m, Angle=Vec2(0.0, 0.0)*NoUnits, W=2000.0m, H=100.0m)
d3 = hostrock_displacement(c3, Point3(200.0, 200.0, -5010.0))
@test d3[1] ≈ 0.0
@test d3[2] ≈ 0.0
@test d3[3] ≈ -100.0

d3o = hostrock_displacement(c3, Point3(1200.0, 1200.0, -5010.0))
@test d3o[3] ≈ 0.0

# ---- CylindricalDikeTopAccretionFullModelAdvection ------------------------

cf2 = CylindricalDikeTopAccretionFullModelAdvection(Center=Point2(0.0, -5000.0)*m, Angle=Vec1(0.0)*NoUnits, W=2000.0m, H=100.0m)
@test isdimensional(cf2) == true
@test isdimensional(nondimensionalize(cf2, CharDim)) == false

cf2u = CylindricalDikeTopAccretionFullModelAdvection(cf2, H=80.0m)
@test UnitValue(cf2u.H) ≈ 80.0m

@test InjectSills.area(cf2) ≈ UnitValue(cf2.W) * UnitValue(cf2.H)
@test InjectSills.volume(cf2) ≈ π * (UnitValue(cf2.W) / 2)^2 * UnitValue(cf2.H)

# Full model advection variant: all lower half-space moves down
df_below_anyx = hostrock_displacement(cf2, Point2(50000.0, -5010.0))
@test df_below_anyx[2] ≈ -100.0

df_above = hostrock_displacement(cf2, Point2(0.0, -4990.0))
@test df_above[2] ≈ 0.0

cf3 = CylindricalDikeTopAccretionFullModelAdvection(Center=Point3(0.0, 0.0, -5000.0)*m, Angle=Vec2(0.0, 0.0)*NoUnits, W=2000.0m, H=100.0m)
df3 = hostrock_displacement(cf3, Point3(100.0, 100.0, -5010.0))
@test df3[3] ≈ -100.0

df3o = hostrock_displacement(cf3, Point3(2000.0, 2000.0, -5010.0))
@test df3o[3] ≈ 0.0
