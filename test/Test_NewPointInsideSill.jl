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
        W      = 1000.0m,
        H      = 100.0m,
    ),
    PennyShapedSill(           # rotated
        Center = Point2(0.0, -5000.0) * m,
        Angle  = Vec1(30.0) * NoUnits,
        W      = 1000.0m,
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
        W      = 1000.0m,
        H      = 100.0m,
    ),
    PennyShapedSill(           # rotated dip + strike
        Center = Point3(0.0, 0.0, -5000.0) * m,
        Angle  = Vec2(30.0, 45.0) * NoUnits,
        W      = 1000.0m,
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
