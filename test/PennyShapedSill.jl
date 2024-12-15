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
#@test d[1] ≈ -0.00017132038852454352
#@test d[2] ≈ 0.01078386907125663


# test displacement routines itself 
r = 10
z = 50
ΔP = 1.0
ν = 0.25
E = 1.0e11
W = 0.1


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


#=
using Plots
heatmap(x/1e3,z/1e3,Ux)
=#