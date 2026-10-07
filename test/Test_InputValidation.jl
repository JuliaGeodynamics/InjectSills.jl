using Test
using InjectSills

C2 = Point2(0.0, -5000.0)*m
C3 = Point3(0.0, 0.0, -5000.0)*m

@testset "non-positive dimensions" begin
    for S in (SquareDike, SquareDikeTopAccretion, CylindricalDikeTopAccretion,
              CylindricalDikeTopAccretionFullModelAdvection, EllipticalIntrusion)
        @test_throws "`W` must be positive" S(Center=C2, Angle=Vec1(0.0)*NoUnits, W=-1.0m, H=100.0m)
        @test_throws "`H` must be positive" S(Center=C2, Angle=Vec1(0.0)*NoUnits, W=2000.0m, H=0.0m)
    end
    @test_throws "`R` must be positive" PennyShapedSill(Center=C2, Angle=Vec1(0.0), R=-1000.0m, H=10.0m)
    @test_throws "`r` must be positive" MogiSphere(Center=C2, r=-1500.0m)
    @test_throws "`r` must be positive" McTigueSphere(Center=C2, r=0.0m)
    @test_throws "`ay` must be positive" FiniteEllipsoidalCavity(Center=C3, ax=500.0m, ay=-500.0m, az=2000.0m)
end

@testset "Poisson's ratio" begin
    @test_throws "Poisson's ratio must satisfy -1 < ν < 0.5" PennyShapedSill(Center=C2, Angle=Vec1(0.0), R=1000.0m, H=10.0m, ν=0.5*NoUnits)
    @test_throws "Poisson's ratio must satisfy -1 < ν ≤ 0.5" MogiSphere(Center=C2, ν=0.6*NoUnits)
    @test_throws "Poisson's ratio must satisfy -1 < ν ≤ 0.5" McTigueSphere(Center=C2, ν=-1.0*NoUnits)
    @test MogiSphere(Center=C2, ν=0.5*NoUnits) isa MogiSphere
end

@testset "over-determined penny" begin
    @test_throws "specify at most two of R, H, ΔP, Q" PennyShapedSill(Center=C2, Angle=Vec1(0.0), R=1000.0m, Q=1e6m^3, ΔP=1e6Pa)
    @test_throws "specify at most two of R, H, ΔP, Q" PennyShapedSill(Center=C2, Angle=Vec1(0.0), R=1000.0m)
end

@testset "penny and plane-strain sills take R, not W" begin
    p = PennyShapedSill(Center=C2, Angle=Vec1(0.0), R=1000.0m, H=10.0m)
    s = PlaneStrainSill(Center=C2, Angle=Vec1(0.0), R=1000.0m, H=10.0m)
    @test_throws MethodError PennyShapedSill(Center=C2, Angle=Vec1(0.0), W=1000.0m, H=10.0m)
    @test_throws MethodError PennyShapedSill(p; W=500.0m)
    @test_throws MethodError update_abstractsill(p; W=500.0m)
    @test_throws MethodError PlaneStrainSill(Center=C2, Angle=Vec1(0.0), W=1000.0m, H=10.0m)
    @test_throws MethodError PlaneStrainSill(s; W=500.0m)
    @test_throws MethodError update_abstractsill(s; W=500.0m)
end

@testset "update_abstractsill keywords" begin
    sills = (PennyShapedSill(Center=C2, Angle=Vec1(0.0), R=1000.0m, H=10.0m),
             PlaneStrainSill(Center=C2, Angle=Vec1(0.0), R=1000.0m, H=10.0m),
             SquareDike(Center=C2, Angle=Vec1(0.0)*NoUnits),
             SquareDikeTopAccretion(Center=C2, Angle=Vec1(0.0)*NoUnits),
             CylindricalDikeTopAccretion(Center=C2, Angle=Vec1(0.0)*NoUnits),
             CylindricalDikeTopAccretionFullModelAdvection(Center=C2, Angle=Vec1(0.0)*NoUnits),
             EllipticalIntrusion(Center=C2, Angle=Vec1(0.0)*NoUnits),
             MogiSphere(Center=C2), McTigueSphere(Center=C2),
             FiniteEllipsoidalCavity(Center=C3, ax=500.0m, ay=500.0m, az=2000.0m))
    for s in sills
        @test_throws MethodError update_abstractsill(s; Hieght=1.0m)
    end
end
