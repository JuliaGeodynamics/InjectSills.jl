using Test
using GeoParams, InjectSills

CharDim = GEO_units(length=1000m, temperature=1000C, stress=10Pa, viscosity=1e20Pas)

# ---- construction and basic properties ------------------------------------

s2 = SquareDike(Center=Point2(0.0, -5000.0)*m, Angle=Vec1(0.0)*NoUnits, W=2000.0m, H=100.0m)
@test isdimensional(s2) == true
@test UnitValue(s2.Lengthscale) ≈ 100.0m
@test UnitValue(s2.BoundingBox[1]) ≈ Point2(-1000.0m, -5050.0m)
@test UnitValue(s2.BoundingBox[2]) ≈ Point2(1000.0m, -4950.0m)

s2_nd = nondimensionalize(s2, CharDim)
@test isdimensional(s2_nd) == false

# copy-constructor updates
s2_u = SquareDike(s2, W=2500.0m, H=120.0m)
@test UnitValue(s2_u.W) ≈ 2500.0m
@test UnitValue(s2_u.H) ≈ 120.0m

s2_u2 = SquareDike(s2, W=1800.0)
@test UnitValue(s2_u2.W) ≈ 1800.0m

# geometric consistency
@test InjectSills.area(s2) ≈ UnitValue(s2.W) * UnitValue(s2.H)
@test InjectSills.volume(s2) ≈ UnitValue(s2.W)^2 * UnitValue(s2.H)

# ---- displacement ---------------------------------------------------------

d_above = hostrock_displacement(s2, Point2(0.0, -5000.0 + 10.0))
@test d_above[1] ≈ 0.0
@test d_above[2] ≈ 50.0

d_below = hostrock_displacement(s2, Point2(0.0, -5000.0 - 10.0))
@test d_below[1] ≈ 0.0
@test d_below[2] ≈ -50.0

d_out = hostrock_displacement(s2, Point2(2000.0, -5000.0))
@test d_out[1] ≈ 0.0
@test d_out[2] ≈ 0.0

# array form from generic AbstractSill utility
x = range(-1500.0, 1500.0, length=11)
z = fill(-5000.0, length(x))
dx, dz = hostrock_displacement(s2, collect(x), z)
@test all(dx .≈ 0.0)
@test dz[1] ≈ 0.0
@test dz[end] ≈ 0.0

# rotated 2D case should produce x and z components
s2_rot = SquareDike(Center=Point2(0.0, -5000.0)*m, Angle=Vec1(45.0)*NoUnits, W=2000.0m, H=100.0m)
d_rot = hostrock_displacement(s2_rot, Point2(0.0, -4990.0))
@test abs(d_rot[1]) > 1.0
@test abs(d_rot[2]) > 1.0

# rotate keyword: default applies rotation, rotate=false skips it
p_rot_frame = Point2(900 / sqrt(2), -5000.0 - 900 / sqrt(2))
@test inside(p_rot_frame, s2_rot) == true
@test inside(p_rot_frame, s2_rot; rotate=false) == false

# ---- inside ---------------------------------------------------------------

@test inside(Point2(0.0, -5000.0), s2) == true
@test inside(Point2(999.0, -5000.0), s2) == true
@test inside(Point2(1001.0, -5000.0), s2) == false
@test inside(Point2(0.0, -4950.0), s2) == true
@test inside(Point2(0.0, -4949.0), s2) == false

# ---- 3D -------------------------------------------------------------------

s3 = SquareDike(Center=Point3(0.0, 0.0, -5000.0)*m, Angle=Vec2(0.0, 0.0)*NoUnits, W=2000.0m, H=100.0m)
@test UnitValue(s3.Lengthscale) ≈ 100.0m
@test UnitValue(s3.BoundingBox[1]) ≈ Point3(-1000.0m, -1000.0m, -5050.0m)
@test UnitValue(s3.BoundingBox[2]) ≈ Point3(1000.0m, 1000.0m, -4950.0m)

@test inside(Point3(0.0, 0.0, -5000.0), s3) == true
@test inside(Point3(999.0, 999.0, -5000.0), s3) == true
@test inside(Point3(1001.0, 0.0, -5000.0), s3) == false
@test inside(Point3(0.0, 0.0, -4950.0), s3) == true
@test inside(Point3(0.0, 0.0, -4949.0), s3) == false

d3_above = hostrock_displacement(s3, Point3(0.0, 0.0, -4990.0))
@test d3_above[1] ≈ 0.0
@test d3_above[2] ≈ 0.0
@test d3_above[3] ≈ 50.0

d3_below = hostrock_displacement(s3, Point3(0.0, 0.0, -5010.0))
@test d3_below[3] ≈ -50.0

d3_out = hostrock_displacement(s3, Point3(2000.0, 0.0, -5000.0))
@test d3_out[1] ≈ 0.0
@test d3_out[2] ≈ 0.0
@test d3_out[3] ≈ 0.0
