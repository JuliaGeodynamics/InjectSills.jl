using Test
using GeoParams, InjectSills

CharDim = GEO_units(length=1000m, temperature=1000C, stress=10Pa, viscosity=1e20Pas)

# ---- construction, units, geometry ----------------------------------------
s2 = MindlinSphere(Center=Point2(0.0,-5000.0)*m, r=1500.0m, ΔP=10e6Pa, G=10e9Pa, ν=0.25*NoUnits)
@test isdimensional(s2) == true
@test UnitValue(s2.Lengthscale) ≈ 1500.0m
@test UnitValue(s2.BoundingBox[1]) ≈ Point2(-1500.0m, -6500.0m)
@test isdimensional(nondimensionalize(s2, CharDim)) == false
s2u = MindlinSphere(s2, r=1200.0m, ΔP=12e6Pa)
@test UnitValue(s2u.r) ≈ 1200.0m && UnitValue(s2u.ΔP) ≈ 12e6Pa
@test InjectSills.area(s2) ≈ π * UnitValue(s2.r)^2
@test InjectSills.volume(s2) ≈ 4/3 * π * UnitValue(s2.r)^3
@test inside(Point2(100.0, -5000.0), s2) && !inside(Point2(0.0, -2000.0), s2)
@test UnitValue(update_abstractsill(s2; ν = 0.3 * NoUnits).ν) ≈ 0.3
@test_throws Exception MindlinSphere(Center=Point2(0.0,-5000.0)*m, r=-1.0m)

# ---- displacement field -----------------------------------------------------
r = 1500.0; ΔP = 10e6; G = 10e9; ν = 0.25; c = 5000.0
ΔV = π * ΔP * r^3 / G
s3 = MindlinSphere(Center=Point3(0.0, 0.0, -c)*m, r=r*m, ΔP=ΔP*Pa, G=G*Pa, ν=ν*NoUnits)
mogi3 = MogiSphere(Center=Point3(0.0, 0.0, -c)*m, r=r*m, ΔP=ΔP*Pa, G=G*Pa, ν=ν*NoUnits)

# at the free surface it is the Mogi solution
for p in (Point3(3000.0, 0.0, 0.0), Point3(-1000.0, 2500.0, 0.0), Point3(0.0, 0.0, 0.0))
    @test hostrock_displacement(s3, p) ≈ hostrock_displacement(mogi3, p) rtol = 1e-12
end
@test hostrock_displacement(s2, Point2(3000.0, 0.0)) ≈ hostrock_displacement(MogiSphere(Center=Point2(0.0,-c)*m, r=r*m, ΔP=ΔP*Pa, G=G*Pa, ν=ν*NoUnits), Point2(3000.0, 0.0)) rtol = 1e-12

# close to the source it is the full-space field ΔV/(4π) r/|r|³ (outward)
u0 = ΔV / (4π) / 10.0^2
for (p, e) in ((Point3(10.0, 0.0, -c), (1, 0, 0)), (Point3(0.0, 0.0, -c - 10.0), (0, 0, -1)))
    u = hostrock_displacement(s3, p)
    @test u[1] ≈ e[1] * u0 atol = 1e-3 * u0
    @test u[3] ≈ e[3] * u0 atol = 1e-3 * u0
end

# the free surface is traction free: σ_zz = σ_xz = σ_yz = 0 at z = 0 (finite differences, one-sided in z)
let λ = 2G * ν / (1 - 2ν), h = 0.02
    for (x, y) in ((1000.0, 0.0), (3000.0, 1500.0), (-6000.0, 2000.0))
        Gr = zeros(3, 3)
        for j in 1:3
            dp = zeros(3); dp[j] = h
            zp = min(dp[3], 0.0); zm = -dp[3]
            up = hostrock_displacement(s3, Point3(x + dp[1], y + dp[2], zp)); um = hostrock_displacement(s3, Point3(x - dp[1], y - dp[2], zm))
            st = j == 3 ? (zp - zm) : 2h
            for i in 1:3; Gr[i, j] = (up[i] - um[i]) / st; end
        end
        e = 0.5 .* (Gr .+ Gr'); tr = e[1, 1] + e[2, 2] + e[3, 3]
        scale = 2G * maximum(abs.(e))
        @test abs(λ * tr + 2G * e[3, 3]) < 1e-4 * scale
        @test abs(2G * e[1, 3]) < 1e-4 * scale && abs(2G * e[2, 3]) < 1e-4 * scale
    end
end

# the volume displaced through a box around the source is the volume change ΔV (Mogi: 4(1−ν) ΔV)
let L = 1000.0, n = 120, h = 2L / n, flux = 0.0, flux_mogi = 0.0
    for a in range(-L + h/2, L - h/2, length = n), b in range(-L + h/2, L - h/2, length = n)
        for (s, acc) in ((s3, :m), (mogi3, :g))
            f = (hostrock_displacement(s, Point3(L, a, -c + b))[1] - hostrock_displacement(s, Point3(-L, a, -c + b))[1] +
                 hostrock_displacement(s, Point3(a, L, -c + b))[2] - hostrock_displacement(s, Point3(a, -L, -c + b))[2] +
                 hostrock_displacement(s, Point3(a, b, -c + L))[3] - hostrock_displacement(s, Point3(a, b, -c - L))[3]) * h^2
            acc == :m ? (flux += f) : (flux_mogi += f)
        end
    end
    @test flux ≈ ΔV rtol = 0.01
    @test flux_mogi ≈ 4 * (1 - ν) * ΔV rtol = 0.01
end

# array interface (2D and 3D)
X = [x for x in -3000.0:500:3000, z in -8000.0:500:0]; Z = [z for x in -3000.0:500:3000, z in -8000.0:500:0]
dX, dZ = hostrock_displacement(s2, X, Z)
@test size(dX) == size(X) && dZ[7, 17] ≈ hostrock_displacement(s2, Point2(0.0, 0.0))[2]
