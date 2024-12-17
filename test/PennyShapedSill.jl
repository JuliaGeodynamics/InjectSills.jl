using Test

using GeoParams, InjectSills
CharDim = GEO_units(length=1000m, temperature=1000C, stress=10Pa, viscosity=1e20Pas)

a = GeoUnit(Point2(100.0m))
nondimensionalize(a, CharDim)
a_nd = nondimensionalize(a, CharDim)
@test a_nd.val == Point2(0.1)

# Define sill with different parameters
sill0 = PennyShapedSill()

sill = PennyShapedSill(ΔP = 1e6Pa)
@test UnitValue(sill.W) ≈ UnitValue(sill0.W) 

sill = PennyShapedSill(ΔP = 1e6Pa, Q=1000m^3)
@test sill.H.val ≈ sill0.H.val 

sill = PennyShapedSill(W = UnitValue(sill0.W), H = UnitValue(sill0.H))
@test UnitValue(sill.ΔP) ≈ UnitValue(sill0.ΔP) 

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


# test displacement routines itself 

# test in 2D
sill2D = PennyShapedSill(Center=Point2(0,-15000)*m, H=100.0m, W=10000.0m, Angle=Vec1(80))

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
sill3D = PennyShapedSill(Center=Point3(0.0,0,-25000)*m, H=100.0m, W=20000.0m, Angle=Vec2(0,0))

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
@test all(extrema(Uy) .≈ (-11.219951034388439, 10.885274888931045))
@test all(extrema(Uz) .≈ (-49.09081409472374, 49.99999999998878))

# Perform computations for a few selected points, which we set in MTK by hand
p = Point2(0,12.0)
sill2D = PennyShapedSill(Center=Point2(0,-15000)*m, H=100.0m, W=10000.0m, Angle=Vec1(0))
d = hostrock_displacement(sill2D, p)
@test d[1] ≈ 4.162983686155986e-12 # compared with MTK
@test d[2] ≈ 12.660339641108203

p = Point2(0,12.0)
sill2D = PennyShapedSill(Center=Point2(0,-15000)*m, H=100.0m, W=10000.0m, Angle=Vec1(80))
d = hostrock_displacement(sill2D, p)
@test d[1] ≈ 1.1159342406236077  # compared with MTK
@test d[2] ≈ -1.179429186779476

sill2D = PennyShapedSill(Center=Point2(0,-25000)*m, H=100.0m, W=10000.0m, Angle=Vec1(80))
p = Point2(0,12.0)
d = hostrock_displacement(sill2D, p)
@test d[1] ≈ 0.29345422900633555  # compared with MTK
@test d[2] ≈ -0.5238221250191283

sill3D  = PennyShapedSill(Center=Point3(0.0,0,-25000)*m, H=100.0m, W=10000.0m, Angle=Vec2(0,0))
p       = Point3(0,0.0,12)
d       = hostrock_displacement(sill3D, p)
@test d[1] ≈ 0.0  # compared with MTK
@test d[2] ≈ 0.0
@test d[3] ≈ 5.617623561481996

sill3D  = PennyShapedSill(Center=Point3(0.0,0,-25000)*m, H=100.0m, W=10000.0m, Angle=Vec2(80,0))
p       = Point3(0,0.0,12)
d       = hostrock_displacement(sill3D, p)
@test d[1] ≈ 0.29345422900633555  # compared with MTK
@test d[2] ≈ 0.0
@test d[3] ≈ -0.5238221250191283

sill3D  = PennyShapedSill(Center=Point3(0.0,0,-25000)*m, H=100.0m, W=10000.0m, Angle=Vec2(80.0,-31))
p       = Point3(0,0,12.0)
d       = hostrock_displacement(sill3D, p)
@test d[1] ≈ 0.279395551087826  
@test d[2] ≈ 0.291934850929267
@test d[3] ≈-0.4440914005320243

# Case that used to give NaN:
sill3D  = PennyShapedSill(Center=Point3(0.0,0,-25000)*m, H=100.0m, W=20000.0m, Angle=Vec2(0,0))
p       = Point3(0,-20e3,-25e3)
d       = hostrock_displacement(sill3D, p)
@test d[1] ≈ 0.0
@test d[2] ≈ -11.219951034388439
@test d[3] ≈ 2.2728421032456027e-5

#=
using Plots
heatmap(x/1e3,z/1e3,Ux)
=#