using Test, InteractiveUtils
using GeoParams, InjectSills

# ---------------------------------------------------------------------------
# Generic test for new_point_inside_sill
#
# Strategy: build one representative instance of every concrete AbstractSill
# subtype (2D and 3D), call new_point_inside_sill a handful of times, and
# verify that every returned point is inside the sill via `inside()`.
#
# To cover a new subtype in the future, add an instance to the appropriate
# list below.  The safety-net testset at the bottom will catch any subtype
# that was added to the package but forgotten here.
# ---------------------------------------------------------------------------


# run N_TRIALS random draws per sill to reduce the chance of a false-positive
const N_TRIALS = 100

# ---- 2D sills --------------------------------------------------------------

sills_2d = [
    PennyShapedSill(
        Center = Point2(0.0, -5000.0) * m,
        Angle  = Vec1(0.0) * NoUnits,
        R      = 1000.0m,
        H      = 100.0m,
    ),
    PennyShapedSill(           # rotated
        Center = Point2(0.0, -5000.0) * m,
        Angle  = Vec1(30.0) * NoUnits,
        R      = 1000.0m,
        H      = 100.0m,
    ),
    PlaneStrainSill(           # rotated
        Center = Point2(0.0, -5000.0) * m,
        Angle  = Vec1(30.0) * NoUnits,
        R      = 1000.0m,
        H      = 100.0m,
    ),
    SquareDike(
        Center = Point2(0.0, -5000.0) * m,
        Angle  = Vec1(0.0) * NoUnits,
        W      = 2000.0m,
        H      = 100.0m,
    ),
    SquareDike(                # rotated
        Center = Point2(0.0, -5000.0) * m,
        Angle  = Vec1(45.0) * NoUnits,
        W      = 2000.0m,
        H      = 100.0m,
    ),
    SquareDikeTopAccretion(
        Center = Point2(0.0, -5000.0) * m,
        Angle  = Vec1(0.0) * NoUnits,
        W      = 2000.0m,
        H      = 100.0m,
    ),
    CylindricalDikeTopAccretion(
        Center = Point2(0.0, -5000.0) * m,
        Angle  = Vec1(0.0) * NoUnits,
        W      = 2000.0m,
        H      = 100.0m,
    ),
    CylindricalDikeTopAccretionFullModelAdvection(
        Center = Point2(0.0, -5000.0) * m,
        Angle  = Vec1(0.0) * NoUnits,
        W      = 2000.0m,
        H      = 100.0m,
    ),
    EllipticalIntrusion(
        Center = Point2(0.0, -5000.0) * m,
        Angle  = Vec1(0.0) * NoUnits,
        W      = 2000.0m,
        H      = 200.0m,
    ),
    EllipticalIntrusion(       # rotated
        Center = Point2(0.0, -5000.0) * m,
        Angle  = Vec1(45.0) * NoUnits,
        W      = 2000.0m,
        H      = 200.0m,
    ),
    MogiSphere(
        Center = Point2(0.0, -5000.0) * m,
        r      = 1500.0m,
        ΔP     = 10e6Pa,
        G      = 10e9Pa,
        ν      = 0.25 * NoUnits,
    ),
    McTigueSphere(
        Center = Point2(0.0, -5000.0) * m,
        r      = 1500.0m,
        ΔP     = 10e6Pa,
        G      = 10e9Pa,
        ν      = 0.25 * NoUnits,
    ),
]

@testset "new_point_inside_sill – 2D" begin
    for sill in sills_2d
        @testset "$(nameof(typeof(sill)))" begin
            for _ in 1:N_TRIALS
                pt = new_point_inside_sill(sill)
                @test inside(pt, sill) == true
            end
        end
    end
end

# ---- 3D sills --------------------------------------------------------------

sills_3d = [
    PennyShapedSill(
        Center = Point3(0.0, 0.0, -5000.0) * m,
        Angle  = Vec2(0.0, 0.0) * NoUnits,
        R      = 1000.0m,
        H      = 100.0m,
    ),
    PennyShapedSill(           # rotated dip + strike
        Center = Point3(0.0, 0.0, -5000.0) * m,
        Angle  = Vec2(30.0, 45.0) * NoUnits,
        R      = 1000.0m,
        H      = 100.0m,
    ),
    SquareDike(
        Center = Point3(0.0, 0.0, -5000.0) * m,
        Angle  = Vec2(0.0, 0.0) * NoUnits,
        W      = 2000.0m,
        H      = 100.0m,
    ),
    SquareDike(                # rotated
        Center = Point3(0.0, 0.0, -5000.0) * m,
        Angle  = Vec2(30.0, 0.0) * NoUnits,
        W      = 2000.0m,
        H      = 100.0m,
    ),
    SquareDikeTopAccretion(
        Center = Point3(0.0, 0.0, -5000.0) * m,
        Angle  = Vec2(0.0, 0.0) * NoUnits,
        W      = 2000.0m,
        H      = 100.0m,
    ),
    CylindricalDikeTopAccretion(
        Center = Point3(0.0, 0.0, -5000.0) * m,
        Angle  = Vec2(0.0, 0.0) * NoUnits,
        W      = 2000.0m,
        H      = 100.0m,
    ),
    CylindricalDikeTopAccretionFullModelAdvection(
        Center = Point3(0.0, 0.0, -5000.0) * m,
        Angle  = Vec2(0.0, 0.0) * NoUnits,
        W      = 2000.0m,
        H      = 100.0m,
    ),
    EllipticalIntrusion(
        Center = Point3(0.0, 0.0, -5000.0) * m,
        Angle  = Vec2(0.0, 0.0) * NoUnits,
        W      = 2000.0m,
        H      = 200.0m,
    ),
    EllipticalIntrusion(       # rotated
        Center = Point3(0.0, 0.0, -5000.0) * m,
        Angle  = Vec2(30.0, 45.0) * NoUnits,
        W      = 2000.0m,
        H      = 200.0m,
    ),
    MogiSphere(
        Center = Point3(0.0, 0.0, -5000.0) * m,
        r      = 1500.0m,
        ΔP     = 10e6Pa,
        G      = 10e9Pa,
        ν      = 0.25 * NoUnits,
    ),
    McTigueSphere(
        Center = Point3(0.0, 0.0, -5000.0) * m,
        r      = 1500.0m,
        ΔP     = 10e6Pa,
        G      = 10e9Pa,
        ν      = 0.25 * NoUnits,
    ),
    FiniteEllipsoidalCavity(
        Center = Point3(0.0, 0.0, -5000.0) * m,
        ax     = 500.0m,
        ay     = 500.0m,
        az     = 2000.0m,
        Angle  = Vec{3}(0.0, 0.0, 0.0) * NoUnits,
        ΔP     = 1e6Pa,
        mu     = 10e9Pa,
        lambda = 10e9Pa,
    ),
    FiniteEllipsoidalCavity(   # rotated
        Center = Point3(0.0, 0.0, -5000.0) * m,
        ax     = 500.0m,
        ay     = 300.0m,
        az     = 2000.0m,
        Angle  = Vec{3}(0.0, 0.0, 45.0) * NoUnits,
        ΔP     = 1e6Pa,
        mu     = 10e9Pa,
        lambda = 10e9Pa,
    ),
]

@testset "new_point_inside_sill – 3D" begin
    for sill in sills_3d
        @testset "$(nameof(typeof(sill)))" begin
            for _ in 1:N_TRIALS
                pt = new_point_inside_sill(sill)
                @test inside(pt, sill) == true
            end
        end
    end
end

# Rotated sills: samples must cover the whole rotated body, not only its unrotated box
@testset "new_point_inside_sill – rotated" begin
    penny = PennyShapedSill(Center=Point3(0.0, 0.0, -5000.0)*m, R=1000.0m, H=1.0m, Angle=Vec2(45.0, 0.0))
    pts   = [new_point_inside_sill(penny) for _ in 1:500]
    @test all(p -> inside(p, penny), pts)
    @test extrema(p[3] for p in pts)[2] - extrema(p[3] for p in pts)[1] > 1000   # spans ±707 m in depth

    fec = FiniteEllipsoidalCavity(Center=Point3(0.0, 0.0, -5000.0)*m, ax=300.0m, ay=500.0m, az=1500.0m,
        Angle=Vec{3}(20.0, 40.0, -30.0)*NoUnits, ΔP=10e6Pa, mu=10e9Pa, lambda=10e9Pa)
    @test all(p -> inside(p, fec), [new_point_inside_sill(fec) for _ in 1:500])
end

# Sills are passed to GPU kernels by value; FiniteEllipsoidalCavity holds arrays, which
# `hostrock_displacement!` moves to the device with Adapt.
@testset "isbits" begin
    for s in vcat(sills_2d, sills_3d)
        s isa FiniteEllipsoidalCavity && continue
        @test isbitstype(typeof(s))
    end
end

# ---- safety net: fail if a new AbstractSill subtype has no test instance ---
@testset "coverage – all AbstractSill subtypes are tested" begin
    registered = Set(nameof(S) for S in subtypes(InjectSills.AbstractSill))
    tested     = Set(nameof(typeof(s)) for s in Iterators.flatten((sills_2d, sills_3d)))
    untested   = setdiff(registered, tested)
    if !isempty(untested)
        @warn "AbstractSill subtypes with no new_point_inside_sill test instance" untested
    end
    @test isempty(untested)
end

# a sill that contains no point of its bounding box
struct EmptySill{B} <: AbstractSill{2, Float64}
    BoundingBox::B
end
InjectSills.inside(::Point2, ::EmptySill; rotate::Bool=true) = false

@testset "new_point_inside_sill – no point found" begin
    @test_throws "no point inside EmptySill after 1000 tries" new_point_inside_sill(EmptySill(PennyShapedSill().BoundingBox))
end
