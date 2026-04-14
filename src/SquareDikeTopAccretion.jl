using Adapt
import Base.show
import GeoParams: isdimensional

export SquareDikeTopAccretion

"""
    SquareDikeTopAccretion{N,_T}

Kinematic top-accretion square dike source in 2D or 3D.
Only the lower half-space is displaced vertically (under-accretion).
"""
struct SquareDikeTopAccretion{N, _T, N1, N2, U1, U2, U3} <: AbstractSill{N, _T}
    Center::GeoUnit{Point{N, _T}, U1}
    Angle::GeoUnit{Vec{N1, _T}, U2}
    W::GeoUnit{_T, U1}
    H::GeoUnit{_T, U1}
    Lengthscale::GeoUnit{_T, U1}
    BoundingBox::Tuple
    RotMat::GeoUnit{SMatrix{N, N, _T, N2}, U3}
end
Adapt.@adapt_structure SquareDikeTopAccretion

isdimensional(s::SquareDikeTopAccretion) = isdimensional(s.W)

"""
    SquareDikeTopAccretion(; Center=Point2(0.0, -5000.0)*m, Angle=Vec1(0.0)*NoUnits, W=2000.0m, H=100.0m)
"""
function SquareDikeTopAccretion(;
    Center = Point2(0.0, -5000.0) * m,
    Angle  = Vec1(0.0) * NoUnits,
    W      = 2000.0m,
    H      = 100.0m,
)
    @assert length(Center) == length(Angle) + 1
    RotMat = RotationMatrix(ustrip.(Angle))
    Cg = convert(GeoUnit, Center)
    Wg = convert(GeoUnit, W)
    Hg = convert(GeoUnit, H)
    Lengthscale = Hg
    BoundingBox = if length(Center) == 2
        unrotated_bounding_box(Cg, Wg.val / 2, Hg.val / 2)
    else
        unrotated_bounding_box(Cg, Wg.val / 2, Wg.val / 2, Hg.val / 2)
    end
    return SquareDikeTopAccretion(Cg, convert(GeoUnit, Angle), Wg, Hg, Lengthscale, BoundingBox, convert(GeoUnit, RotMat))
end

"""
    SquareDikeTopAccretion(s::SquareDikeTopAccretion; kwargs...)
"""
function SquareDikeTopAccretion(s::SquareDikeTopAccretion; kwargs...)
    valid = (:Center, :Angle, :W, :H)
    all(k -> k in valid, keys(kwargs)) ||
        error("Invalid keyword for SquareDikeTopAccretion(s; ...). Valid keys are: $(valid)")

    base = (
        Center = UnitValue(s.Center),
        Angle  = UnitValue(s.Angle),
        W      = UnitValue(s.W),
        H      = UnitValue(s.H),
    )

    kw = Dict{Symbol, Any}(kwargs)
    for sym in (:W, :H)
        if haskey(kw, sym) && kw[sym] isa Number && !(kw[sym] isa typeof(oneunit(getproperty(base, sym))))
            kw[sym] = kw[sym] * oneunit(getproperty(base, sym))
        end
    end

    return SquareDikeTopAccretion(; merge(base, (; kw...))...)
end

function show(io::IO, s::SquareDikeTopAccretion)
    label = isdimensional(s) ? "dimensional units" : "nondimensional"
    println(io, "Square dike top-accretion ($label):")
    println(io, "   Center      : $((s.Center.val...,).*s.Center.unit)")
    println(io, "   Angle       : $(s.Angle.val)")
    println(io, "   Width (W)   : $(UnitValue(s.W))")
    println(io, "   Thickness(H): $(UnitValue(s.H))")
    return nothing
end

area(s::SquareDikeTopAccretion) = UnitValue(s.W) * UnitValue(s.H)
volume(s::SquareDikeTopAccretion) = UnitValue(s.W)^2 * UnitValue(s.H)

"""
    d = hostrock_displacement(sill::SquareDikeTopAccretion{N,_T}, p::Point{N,_T})

Kinematic opening field used by MTK `Type="SquareDike_TopAccretion"`:
- inside footprint and below center plane: vertical displacement `-H`
- elsewhere: zero displacement
"""
function hostrock_displacement(sill::SquareDikeTopAccretion{N, _T}, p::Point{N, _T}) where {N, _T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill

    p_r = rotate_point(p - Center, RotMat)

    if N == 2
        if abs(p_r[1]) <= W / 2 && p_r[2] <= 0
            d_r = Vec2{_T}(zero(_T), -H)
            return rotate_point(d_r, RotMat')
        end
        return Vec2{_T}(zero(_T), zero(_T))
    else
        if abs(p_r[1]) <= W / 2 && abs(p_r[2]) <= W / 2 && p_r[3] <= 0
            d_r = Vec3{_T}(zero(_T), zero(_T), -H)
            return rotate_point(d_r, RotMat')
        end
        return Vec3{_T}(zero(_T), zero(_T), zero(_T))
    end
end

function inside(p::Point{2, _T}, sill::SquareDikeTopAccretion{2, _T}; rotate::Bool=true) where {_T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill
    p_r = p - Center
    if rotate
        p_r = rotate_point(p_r, RotMat)
    end
    return abs(p_r[1]) <= W / 2 && abs(p_r[2]) <= H / 2
end

function inside(p::Point{3, _T}, sill::SquareDikeTopAccretion{3, _T}; rotate::Bool=true) where {_T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill
    p_r = p - Center
    if rotate
        p_r = rotate_point(p_r, RotMat)
    end
    return abs(p_r[1]) <= W / 2 && abs(p_r[2]) <= W / 2 && abs(p_r[3]) <= H / 2
end

function update_abstractsill(s::SquareDikeTopAccretion; kwargs...)
    params = (
        Center = UnitValue(s.Center),
        Angle  = UnitValue(s.Angle),
        W      = UnitValue(s.W),
        H      = UnitValue(s.H),
    )
    return SquareDikeTopAccretion(; merge(params, kwargs)...)
end
