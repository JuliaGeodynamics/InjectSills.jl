using Test
using InjectSills

sill2 = PennyShapedSill(Center=Point2(0.0, -5000.0)*m, Angle=Vec1(30.0), R=1000.0m, H=10.0m)
sill3 = PennyShapedSill(Center=Point3(0.0, 0.0, -5000.0)*m, Angle=Vec2(20.0, 40.0), R=1000.0m, H=10.0m)
X = [x for x in range(-2000.0, 2000.0, length=7), _ in 1:5]
Y = [y for _ in 1:7, y in range(-1500.0, 1500.0, length=5)]
Z = fill(-4900.0, 7, 5)

@testset "matches the point method" begin
    Dx, Dz = hostrock_displacement(sill2, X, Z)
    @test all(I -> (Dx[I], Dz[I]) == Tuple(hostrock_displacement(sill2, Point2(X[I], Z[I]))), eachindex(X))
    Dx, Dy, Dz = hostrock_displacement(sill3, X, Y, Z)
    @test all(I -> (Dx[I], Dy[I], Dz[I]) == Tuple(hostrock_displacement(sill3, Point3(X[I], Y[I], Z[I]))), eachindex(X))
end

@testset "FiniteEllipsoidalCavity" begin
    fec = FiniteEllipsoidalCavity(Center=Point3(0.0, 0.0, -5000.0)*m, ax=500.0m, ay=300.0m, az=2000.0m,
        Angle=Vec{3}(10.0, 20.0, 45.0)*NoUnits, ΔP=1e6Pa, mu=10e9Pa, lambda=10e9Pa)
    Ux, Uy, Uz = hostrock_displacement(fec, X, Y, zero(X))
    @test (Ux, Uy, Uz) == hostrock_displacement(fec, X, Y)[1:3]
    @test_throws "all z coordinates must be 0" hostrock_displacement(fec, X, Y, fill(-10.0, size(X)))
end

@testset "views" begin
    ref = hostrock_displacement(sill3, X, Y, Z)
    Dv  = hostrock_displacement(sill3, view(X, 2:6, :), view(Y, 2:6, :), view(Z, 2:6, :))
    @test all(k -> Dv[k] == ref[k][2:6, :], 1:3)
end

@testset "hostrock_displacement!" begin
    @test_throws "all arrays must have axes" hostrock_displacement!((similar(X), similar(X)), sill2, (X, Z[1:6, :]))
    Xn = copy(X); Xn[3] = NaN
    Dx, Dz = hostrock_displacement!((similar(X), similar(X)), sill2, (Xn, Z); skipnan=true)
    @test Dx[3] == Dz[3] == 0
    @test Dx[4] == hostrock_displacement(sill2, Point2(X[4], Z[4]))[1]
end

@testset "Float32" begin
    function sills(::Type{F}) where {F}
        C2, C3 = Point2{F}(0, -5000)*m, Point3{F}(0, 0, -5000)*m
        (PennyShapedSill(Center=C2, Angle=Vec1{F}(30)*NoUnits, R=F(1000)m, H=F(10)m, E=F(1.5e10)Pa, ν=F(0.3)*NoUnits),
         PennyShapedSill(Center=C3, Angle=Vec2{F}(30, 45)*NoUnits, R=F(1000)m, H=F(10)m, E=F(1.5e10)Pa, ν=F(0.3)*NoUnits),
         EllipticalIntrusion(Center=C3, Angle=Vec2{F}(30, 45)*NoUnits, W=F(2000)m, H=F(200)m),
         McTigueSphere(Center=C3, r=F(1500)m, ΔP=F(10e6)Pa, G=F(10e9)Pa, ν=F(0.25)*NoUnits),
         FiniteEllipsoidalCavity(Center=C3, ax=F(500)m, ay=F(300)m, az=F(2000)m, Angle=Vec{3,F}(10, 20, 45)*NoUnits,
                                 ΔP=F(1e6)Pa, mu=F(10e9)Pa, lambda=F(10e9)Pa))
    end
    for (s32, s64) in zip(sills(Float32), sills(Float64))
        N  = length(s32.Center.val)
        x32, x64 = Float32.(collect(range(-3000, 3000, length=6))), collect(range(-3000.0, 3000.0, length=6))
        z = s32 isa FiniteEllipsoidalCavity ? 0.0 : -4800.0
        args32 = N == 2 ? (x32, fill(Float32(z), 6)) : (x32, x32 ./ 2, fill(Float32(z), 6))
        args64 = N == 2 ? (x64, fill(z, 6))          : (x64, x64 ./ 2, fill(z, 6))
        D32, D64 = hostrock_displacement(s32, args32...), hostrock_displacement(s64, args64...)
        @test all(D -> eltype(D) == Float32, D32)
        scale = maximum(D -> maximum(abs, D), D64)
        @test all(k -> maximum(abs, D32[k] .- D64[k]) / scale < 1e-5, 1:N)
    end
end
