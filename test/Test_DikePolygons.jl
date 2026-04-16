using Test
using GeoParams, InjectSills

@testset "dike_polygon fallback" begin
    s = SquareDike(W=1000.0m, H=100.0m)
    @test_throws ErrorException dike_polygon(s)
end

@testset "dike_polygon cylindrical top-accretion" begin
    s = CylindricalDikeTopAccretion(Center=Point2(10.0, -5000.0)*m, W=2000.0m, H=100.0m)
    poly = dike_polygon(s, 21)
    @test length(poly) == 2
    x, z = poly
    @test length(x) == 42
    @test length(z) == 42
    @test minimum(x) ≈ -990.0
    @test maximum(x) ≈ 1010.0
    @test first(z) ≈ -4900.0
    @test last(z) ≈ -5050.0
end

@testset "dike_polygon cylindrical top-accretion full-model" begin
    s = CylindricalDikeTopAccretionFullModelAdvection(Center=Point2(0.0, -3000.0)*m, W=1000.0m, H=80.0m)
    poly = dike_polygon(s, 11)
    x, z = poly
    @test length(x) == 22
    @test length(z) == 22
    @test minimum(x) ≈ -500.0
    @test maximum(x) ≈ 500.0
    @test first(z) ≈ -2920.0
    @test last(z) ≈ -3040.0
end

@testset "dike_polygon elliptical intrusion" begin
    s = EllipticalIntrusion(Center=Point2(0.0, -2000.0)*m, W=1000.0m, H=200.0m)
    poly = dike_polygon(s, 101)
    @test length(poly) == 2
    x, z = poly
    @test length(x) == 101
    @test length(z) == 101
    @test maximum(x) ≈ 500.0 atol=1e-8
    @test minimum(x) ≈ -500.0 atol=1e-8
    @test maximum(z) ≈ -1900.0 atol=1e-8
    @test minimum(z) ≈ -2100.0 atol=1e-8
end

@testset "dike_polygon penny-shaped sill" begin
    s = PennyShapedSill(Center=Point2(10.0, -5000.0)*m, W=1000.0m, H=200.0m, Angle=Vec1(0.0)*NoUnits)
    poly = dike_polygon(s, 101)
    x, z = poly
    @test length(x) == 101
    @test length(z) == 101
    @test maximum(x) ≈ 1010.0 atol=1e-8
    @test minimum(x) ≈ -990.0 atol=1e-8
    @test maximum(z) ≈ -4900.0 atol=1e-8
    @test minimum(z) ≈ -5100.0 atol=1e-8

    srot = PennyShapedSill(Center=Point2(10.0, -5000.0)*m, W=1000.0m, H=200.0m, Angle=Vec1(90.0)*NoUnits)
    polyrot = dike_polygon(srot, 101)
    xr, zr = polyrot
    @test maximum(xr) ≈ 110.0 atol=1e-6
    @test minimum(xr) ≈ -90.0 atol=1e-6
    @test maximum(zr) ≈ -4000.0 atol=1e-6
    @test minimum(zr) ≈ -6000.0 atol=1e-6
end

@testset "dike_polygon spheres in 2D" begin
    mogi = MogiSphere(Center=Point2(10.0, -3000.0)*m, r=500.0m)
    x, z = dike_polygon(mogi, 101)
    @test length(x) == 101
    @test length(z) == 101
    @test maximum(x) ≈ 510.0 atol=1e-8
    @test minimum(x) ≈ -490.0 atol=1e-8
    @test maximum(z) ≈ -2500.0 atol=1e-8
    @test minimum(z) ≈ -3500.0 atol=1e-8

    mct = McTigueSphere(Center=Point2(-20.0, -2000.0)*m, r=400.0m)
    x2, z2 = dike_polygon(mct, 101)
    @test length(x2) == 101
    @test length(z2) == 101
    @test maximum(x2) ≈ 380.0 atol=1e-8
    @test minimum(x2) ≈ -420.0 atol=1e-8
    @test maximum(z2) ≈ -1600.0 atol=1e-8
    @test minimum(z2) ≈ -2400.0 atol=1e-8
end
