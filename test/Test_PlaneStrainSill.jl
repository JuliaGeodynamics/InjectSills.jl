using Test, LinearAlgebra, InteractiveUtils
using GeoParams, InjectSills
include("field_checks.jl")

CharDim = GEO_units(length=1000m, temperature=1000C, stress=10Pa, viscosity=1e20Pas)

# Plane-strain pressurized crack on y = 0, |x| < a, with p/(2G) = H/(4(1-ν)a), in the
# dimensional Westergaard form; returns (u_x, u_y) for y > 0.
function crack_ref(x, y, a, H, ν)
    ζ = complex(x, y); s = sqrt(ζ - a) * sqrt(ζ + a)
    Ẑ = -a^2 / (ζ + s); Z = a^2 / (s * (ζ + s)); k = H / (4 * (1 - ν) * a)
    return k * ((1 - 2ν) * real(Ẑ) - y * imag(Z)), k * (2 * (1 - ν) * imag(Ẑ) - y * real(Z))
end

a, ν, E, ΔP = 2000.0, 0.25, 1.5e10, 3e6

@testset "construction" begin
    s = PlaneStrainSill(Center=Point2(0.0, -5000.0)*m, R=a*m, ΔP=ΔP*Pa, E=E*Pa, ν=ν*NoUnits)
    H = UnitValue(s.H)
    @test H ≈ 4 * (1 - ν^2) * ΔP * a / E * m
    @test UnitValue(s.Q) ≈ π * a * m * H / 2
    @test InjectSills.area(s) ≈ UnitValue(s.Q)
    @test InjectSills.volume(s) ≈ InjectSills.volume(PennyShapedSill(R=a*m, H=H))
    # every pair gives the same sill
    for kw in ((R=a*m, H=H), (R=a*m, Q=s.Q.val*m^2), (H=H, Q=s.Q.val*m^2), (H=H, ΔP=ΔP*Pa), (Q=s.Q.val*m^2, ΔP=ΔP*Pa))
        t = PlaneStrainSill(; Center=Point2(0.0, -5000.0)*m, E=E*Pa, ν=ν*NoUnits, kw...)
        @test all(f -> getfield(t, f).val ≈ getfield(s, f).val, (:R, :H, :ΔP, :Q))
    end
    @test UnitValue(PlaneStrainSill().Q) ≈ 1000.0m^2
    @test UnitValue(PlaneStrainSill().ΔP) ≈ 1e6Pa
    @test UnitValue(PlaneStrainSill(ΔP=2e6Pa).R) ≈ UnitValue(PlaneStrainSill().R)
    @test UnitValue(update_abstractsill(s, R=1000.0m).H) == H
    @test UnitValue(PlaneStrainSill(s, ΔP=6e6).R) ≈ a / 2 * m
    @test isbitstype(typeof(s))
    @test nondimensionalize(s, CharDim) isa PlaneStrainSill
    @test occursin("Plane-strain sill", sprint(show, s))
    @test length(dike_polygon(s)[1]) == 101
    @test_throws "`Center` must be 2D" PlaneStrainSill(Center=Point3(0.0, 0.0, 0.0)*m, R=a*m, H=1.0m)
    @test_throws "specify at most two of R, H, ΔP, Q" PlaneStrainSill(R=a*m, H=1.0m, ΔP=ΔP*Pa)
    @test_throws "`H` must be positive" PlaneStrainSill(R=a*m, H=-1.0m)
    @test_throws MethodError update_abstractsill(s; Hieght=1.0m)
end

@testset "displacement" begin
    s = PlaneStrainSill(Center=Point2(0.0, 0.0)*m, R=a*m, ΔP=ΔP*Pa, E=E*Pa, ν=ν*NoUnits)
    H = s.H.val
    # opening (H/2)√(1 - x²/a²) per face; the point y = 0 is the upper face
    for x in (0.0, 500.0, 1999.0)
        @test hostrock_displacement(s, Point2(x, 0.0))[2] ≈ H / 2 * sqrt(1 - x^2 / a^2)
        @test hostrock_displacement(s, Point2(x, -3.0))[2] == -hostrock_displacement(s, Point2(x, 3.0))[2]
    end
    @test all(isfinite, hostrock_displacement(s, Point2(a, 0.0)))     # tip
    # the faces carry the traction -ΔP
    u(x) = collect(hostrock_displacement(s, Point2(x...)))
    G = E / (2 * (1 + ν))
    for x in (0.0, 1500.0)
        σ = G * stress(u, [x, 1e-3], 1e-4, ν)
        @test σ[2, 2] ≈ -ΔP rtol=1e-6
        @test abs(σ[1, 2]) < 1e-6ΔP
    end
    # area/2 crosses every horizontal line above the sill, -area/2 every line below
    for h in (1.0, 500.0, 5000.0, -1.0, -3000.0)
        @test line_integral(x -> hostrock_displacement(s, Point2(x, h))[2], 20a; breaks=(a,)) ≈ sign(h) * InjectSills.area(s).val / 2 rtol=1e-8
    end
    # rotated sill: the Westergaard field in the sill frame
    sr = PlaneStrainSill(Center=Point2(300.0, -5000.0)*m, Angle=Vec1(25.0), R=a*m, H=H*m, ν=ν*NoUnits)
    R  = sr.RotMat.val
    for q in ([10.0, 11.0], [-2500.0, 800.0], [1000.0, -40.0])
        ux, uy = crack_ref(q[1], abs(q[2]), a, H, ν)
        @test hostrock_displacement(sr, sr.Center.val + R' * Vec2(q...)) ≈ R' * [ux, sign(q[2]) * uy]
    end
    # Float32
    s32 = PlaneStrainSill(Center=Point2f(0, 0)*m, Angle=Vec1f(0)*NoUnits, R=Float32(a)*m, H=Float32(H)*m, E=1.5f10Pa, ν=Float32(ν)*NoUnits)
    @test hostrock_displacement(s32, Point2f(700, 300)) ≈ hostrock_displacement(s, Point2(700.0, 300.0)) rtol=1e-5
    @test !occursin("Float64", string(code_typed(hostrock_displacement, (typeof(s32), Point2f))))
end
