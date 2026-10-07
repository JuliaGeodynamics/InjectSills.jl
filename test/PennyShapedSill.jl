using Test

using GeoParams, InjectSills
CharDim = GEO_units(length=1000m, temperature=1000C, stress=10Pa, viscosity=1e20Pas)

a = GeoUnit(Point2(100.0m))
nondimensionalize(a, CharDim)
a_nd = nondimensionalize(a, CharDim)
@test a_nd.val == Point2(0.1)

# Define sill with different parameters
sill0 = PennyShapedSill()
@test isdimensional(sill0) == true
@test UnitValue(sill0.Lengthscale) ≈ UnitValue(sill0.H)
@test UnitValue(sill0.BoundingBox[1])[1] ≈ UnitValue(sill0.Center)[1] - UnitValue(sill0.R)
@test UnitValue(sill0.BoundingBox[2])[2] ≈ UnitValue(sill0.Center)[2] + UnitValue(sill0.H) / 2

# test that we can nondimensionalize this
sill0_nd = nondimensionalize(sill0, CharDim)
@test isdimensional(sill0_nd) == false

sill = PennyShapedSill(ΔP = 1e6Pa)
@test UnitValue(sill.R) ≈ UnitValue(sill0.R) 

sill_lowP = PennyShapedSill(ΔP = 1e5Pa)
@test UnitValue(sill_lowP.R) ≈ UnitValue(sill0.R)
@test UnitValue(sill_lowP.Q) ≈ 100.0m^3

sill = PennyShapedSill(ΔP = 1e6Pa, Q=1000m^3)
@test sill.H.val ≈ sill0.H.val 

sill = PennyShapedSill(R = UnitValue(sill0.R), H = UnitValue(sill0.H))
@test UnitValue(sill.ΔP) ≈ UnitValue(sill0.ΔP) 

# geometric consistency: planform area A = pi*R^2 and volume Q = 2*pi*H*R^2/3
area0 = π * UnitValue(sill0.R)^2
@test area0 ≈ π * UnitValue(sill0.R)^2
@test UnitValue(sill0.Q) ≈ (2 * π * UnitValue(sill0.H) * UnitValue(sill0.R)^2) / 3

# single-parameter update dispatch from existing sill
sillR = PennyShapedSill(sill0, R = 100.0m)
@test UnitValue(sillR.R) ≈ 100.0m
@test UnitValue(sillR.H) ≈ UnitValue(sill0.H)

sillE = PennyShapedSill(sill0, E = 2.0e10Pa)
@test UnitValue(sillE.E) ≈ 2.0e10Pa

sillRH = PennyShapedSill(sill0, R = 100.0m, H = 50.0m)
@test UnitValue(sillRH.R) ≈ 100.0m
@test UnitValue(sillRH.H) ≈ 50.0m
@test UnitValue(sillRH.Q) ≈ (2 * π * UnitValue(sillRH.H) * UnitValue(sillRH.R)^2) / 3

areaWH = π * UnitValue(sillRH.R)^2
@test areaWH ≈ π * (100.0m)^2

sillMulti = PennyShapedSill(sill0, R = 120.0m, H = 60.0m, E = 2.2e10Pa)
@test UnitValue(sillMulti.R) ≈ 120.0m
@test UnitValue(sillMulti.H) ≈ 60.0m
@test UnitValue(sillMulti.E) ≈ 2.2e10Pa

sillRΔP = PennyShapedSill(sill0, R = 120.0m, ΔP = 2e6Pa)
@test UnitValue(sillRΔP.R) ≈ 120.0m
@test UnitValue(sillRΔP.ΔP) ≈ 2e6Pa
@test UnitValue(sillRΔP.Q) ≈ 16 * (1 - UnitValue(sill0.ν)^2) * UnitValue(sillRΔP.ΔP) * UnitValue(sillRΔP.R)^3 / (3 * UnitValue(sillRΔP.E))

sill2D = PennyShapedSill()
@test length(sill2D.Angle.val) == 1

sill3D = PennyShapedSill(Center=Point3(0.0,0.0,0.0)*m, Angle=Vec2(0.0,0.0))
@test length(sill3D.Angle.val) == 2

sill2D = PennyShapedSill(Center=Point2(0.0,0.0)*m, Angle=Vec1(0.0))
@test length(sill2D.Angle.val) == 1


# displacement:
p = Point2(10.0,11.0)
d = hostrock_displacement(sill0, p)
@test d[1] ≈ -0.00011540383661873777
@test d[2] ≈ 0.010816226990412455

function perform_computation(x::_T,y::_T, sill::AbstractSill{N,_T}) where {N,_T}
    p = Point2{_T}(x,y)
    d = hostrock_displacement(sill, p)
    return d
end
d = perform_computation(10.0, 11.1, sill2D)
@test   -0.00011488707183608771 ≈ d[1]


# test displacement routines itself 

# test in 2D
sill2D = PennyShapedSill(Center=Point2(0,-15000)*m, H=100.0m, R=10000.0m, Angle=Vec1(80))

nx,nz = 129,129
x = range(-15000, stop=15000, length=nx)
z = range(-30000, stop=0, length=nz)
Ux = zeros(nx,nz)
Uz = zeros(nx,nz)

for I in CartesianIndices(Ux)
    Ux[I], Uz[I]  = hostrock_displacement(sill2D, Point2(x[I[1]],z[I[2]]))
end
@test all(extrema(Ux) .≈ (-49.18113369822504, 49.240387455777366))
@test all(extrema(Uz) .≈ (-14.034402454868314, 14.034402454868314))


# test in 3D
sill3D = PennyShapedSill(Center=Point3(0.0,0,-25000)*m, H=100.0m, R=20000.0m, Angle=Vec2(0,0))

nx,ny,nz = 65,65,65
x = range(-15000, stop=15000, length=nx)
y = range(-20000, stop=20000, length=ny)
z = range(-50000, stop=0, length=nz)
Ux = zeros(nx,ny,nz)
Uy = zeros(nx,ny,nz)
Uz = zeros(nx,ny,nz)

for I in CartesianIndices(Ux)
    Ux[I],Uy[I],Uz[I]  = hostrock_displacement(sill3D, Point3(x[I[1]],y[I[2]],z[I[3]]))
end
@test all(extrema(Ux) .≈ (-8.414980322087171, 8.414980322087171))
@test all(extrema(Uy) .≈ (-11.219951034388439, 11.219951034388439))
@test all(extrema(Uz) .≈ (-49.09081409472374, 49.99999999998878))

# Perform computations for a few selected points, which we set in MTK by hand
p = Point2(0,12.0)
sill2D = PennyShapedSill(Center=Point2(0,-15000)*m, H=100.0m, R=10000.0m, Angle=Vec1(0))
d = hostrock_displacement(sill2D, p)
@test d[1] == 0.0                  # on the axis; MTK returns ~4e-12 from its r regularization
@test d[2] ≈ 12.660339641108203

p = Point2(0,12.0)
sill2D = PennyShapedSill(Center=Point2(0,-15000)*m, H=100.0m, R=10000.0m, Angle=Vec1(80))
d = hostrock_displacement(sill2D, p)
@test d[1] ≈ 1.1159342406236077  # compared with MTK
@test d[2] ≈ -1.179429186779476

sill2D = PennyShapedSill(Center=Point2(0,-25000)*m, H=100.0m, R=10000.0m, Angle=Vec1(80))
p = Point2(0,12.0)
d = hostrock_displacement(sill2D, p)
@test d[1] ≈ 0.29345422900633555  # compared with MTK
@test d[2] ≈ -0.5238221250191283

sill3D  = PennyShapedSill(Center=Point3(0.0,0,-25000)*m, H=100.0m, R=10000.0m, Angle=Vec2(0,0))
p       = Point3(0,0.0,12)
d       = hostrock_displacement(sill3D, p)
@test d[1] ≈ 0.0  # compared with MTK
@test d[2] ≈ 0.0
@test d[3] ≈ 5.617623561481996

sill3D  = PennyShapedSill(Center=Point3(0.0,0,-25000)*m, H=100.0m, R=10000.0m, Angle=Vec2(80,0))
p       = Point3(0,0.0,12)
d       = hostrock_displacement(sill3D, p)
@test d[1] ≈ 0.29345422900633555  # compared with MTK
@test d[2] ≈ 0.0
@test d[3] ≈ -0.5238221250191283

sill3D  = PennyShapedSill(Center=Point3(0.0,0,-25000)*m, H=100.0m, R=10000.0m, Angle=Vec2(80.0,-31))
p       = Point3(0,0,12.0)
d       = hostrock_displacement(sill3D, p)
@test d[1] ≈ 0.2515393693569801  
@test d[2] ≈ 0.15114010118163726
@test d[3] ≈ -0.5238221250191283

# Case that used to give NaN:
sill3D  = PennyShapedSill(Center=Point3(0.0,0,-25000)*m, H=100.0m, R=20000.0m, Angle=Vec2(0,0))
p       = Point3(0,-20e3,-25e3)
d       = hostrock_displacement(sill3D, p)
@test d[1] ≈ 0.0
@test d[2] ≈ 11.219951034388439
@test d[3] ≈ 2.2728421032456027e-5

# Axisymmetry in the sill frame: rotating a point 90° about the sill normal rotates its displacement
sill_rot = PennyShapedSill(Center=Point3(0.0,0,-5000)*m, H=10.0m, R=1000.0m, Angle=Vec2(30.0,20.0))
R  = sill_rot.RotMat.val
to_world(q) = Point3(0.0,0,-5000) + R' * q
u1 = R * hostrock_displacement(sill_rot, to_world(Point3(300.0, 400.0, 100.0)))
u2 = R * hostrock_displacement(sill_rot, to_world(Point3(-400.0, 300.0, 100.0)))
@test u2 ≈ [-u1[2], u1[1], u1[3]]

# Case that used to produce a pathological 2D outlier in JustPIC
sill2D_path = PennyShapedSill(Center=Point2(0,-5000)*m, H=40.0m, R=2000.0m, Angle=Vec1(30))
p_path = Point2(153.39749773478215, -4999.988823166007)
d_path = hostrock_displacement(sill2D_path, p_path)
@test isfinite(d_path[1])
@test isfinite(d_path[2])
@test abs(d_path[1]) < 10_000
@test abs(d_path[2]) < 10_000


# test inside routines
sill2D = PennyShapedSill(Center=Point2(0,-15000)*m, H=100.0m, R=10000.0m, Angle=Vec1(0))
p=Point2(0,0.0)
@test inside(p, sill2D) == false

@test inside(Point2(0,-15e3   ),        sill2D) == true
@test inside(Point2(0,-15e3+50),        sill2D) == true
@test inside(Point2(0,-15e3+51),        sill2D) == false
@test inside(Point2(-10000,-15e3    ),  sill2D) == true
@test inside(Point2(-10001,-15e3    ),  sill2D) == false
@test inside(Point2( 10001,-15e3    ),  sill2D) == false


sill3D = PennyShapedSill(Center=Point3(0,0,-15000)*m, H=100.0m, R=10000.0m, Angle=Vec2(0,0))
p=Point3(0,0.0,0)
@test inside(p, sill3D) == false

@test inside(Point3(0,0,-15e3    ),         sill3D) == true
@test inside(Point3(0,0,-15e3+50 ),         sill3D) == true
@test inside(Point3(0,0,-15e3+51 ),         sill3D) == false
@test inside(Point3(-10000,0,-15e3    ),    sill3D) == true
@test inside(Point3(-10001,0,-15e3    ),    sill3D) == false
@test inside(Point3( 10001,0,-15e3    ),    sill3D) == false
@test inside(Point3(0,-10000,-15e3    ),    sill3D) == true
@test inside(Point3(0,-10001,-15e3    ),    sill3D) == false
@test inside(Point3(0, 10001,-15e3    ),    sill3D) == false

#=
using Plots
heatmap(x/1e3,z/1e3,Ux)
=#
