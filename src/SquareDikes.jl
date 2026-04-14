using Adapt
import Base.show
import GeoParams: isdimensional

export SquareDike

"""
    SquareDike{N,_T}

Axis-aligned square/rectangular dike source in 2D or 3D.

Parameters:
====
- `Center::Point{N,_T}` - dike center
- `Angle::Vec{N-1,_T}`  - dip (2D) or dip/strike (3D), in degrees
- `W::_T`               - horizontal width (x and y in 3D)
- `H::_T`               - vertical thickness
- `RotMat`              - precomputed rotation matrix
"""
struct SquareDike{N, _T, N1, N2, U1, U2, U3} <: AbstractSill{N, _T}
    Center::GeoUnit{Point{N, _T}, U1}
    Angle::GeoUnit{Vec{N1, _T}, U2}
    W::GeoUnit{_T, U1}
    H::GeoUnit{_T, U1}
    RotMat::GeoUnit{SMatrix{N, N, _T, N2}, U3}
end
Adapt.@adapt_structure SquareDike

isdimensional(s::SquareDike) = isdimensional(s.W)

"""
    SquareDike(; Center=Point2(0.0, -5000.0)*m, Angle=Vec1(0.0)*NoUnits, W=2000.0m, H=100.0m)

Construct a square-dike source with keyword arguments.
"""
function SquareDike(;
    Center = Point2(0.0, -5000.0) * m,
    Angle  = Vec1(0.0) * NoUnits,
    W      = 2000.0m,
    H      = 100.0m,
)
    @assert length(Center) == length(Angle) + 1
    RotMat = RotationMatrix(ustrip.(Angle))
    return SquareDike(
        convert(GeoUnit, Center),
        convert(GeoUnit, Angle),
        convert(GeoUnit, W),
        convert(GeoUnit, H),
        convert(GeoUnit, RotMat),
    )
end

"""
    SquareDike(s::SquareDike; kwargs...)

Create a new `SquareDike` from an existing one by overriding any subset of
`Center`, `Angle`, `W`, `H`.
"""
function SquareDike(s::SquareDike; kwargs...)
    valid = (:Center, :Angle, :W, :H)
    all(k -> k in valid, keys(kwargs)) ||
        error("Invalid keyword for SquareDike(s; ...). Valid keys are: $(valid)")

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

    return SquareDike(; merge(base, (; kw...))...)
end

function show(io::IO, s::SquareDike)
    label = isdimensional(s) ? "dimensional units" : "nondimensional"
    println(io, "Square dike ($label):")
    println(io, "   Center      : $((s.Center.val...,).*s.Center.unit)")
    println(io, "   Angle       : $(s.Angle.val)")
    println(io, "   Width (W)   : $(UnitValue(s.W))")
    println(io, "   Thickness(H): $(UnitValue(s.H))")
    return nothing
end

area(s::SquareDike) = UnitValue(s.W) * UnitValue(s.H)
volume(s::SquareDike) = UnitValue(s.W)^2 * UnitValue(s.H)

"""
    d = hostrock_displacement(sill::SquareDike{N,_T}, p::Point{N,_T})

Kinematic opening field used by MTK `Type="SquareDike"`:
- inside the square footprint: vertical displacement ±H/2 away from the center plane
- outside: zero displacement
"""
function hostrock_displacement(sill::SquareDike{N, _T}, p::Point{N, _T}) where {N, _T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill

    p_r = rotate_point(p - Center, RotMat)

    if N == 2
        if abs(p_r[1]) <= W / 2
            if p_r[2] < 0
                d_r = Vec2{_T}(zero(_T), -H / 2)
            elseif p_r[2] > 0
                d_r = Vec2{_T}(zero(_T), H / 2)
            else
                d_r = Vec2{_T}(zero(_T), zero(_T))
            end
            return rotate_point(d_r, RotMat')
        end
        return Vec2{_T}(zero(_T), zero(_T))
    else
        if abs(p_r[1]) <= W / 2 && abs(p_r[2]) <= W / 2
            if p_r[3] < 0
                d_r = Vec3{_T}(zero(_T), zero(_T), -H / 2)
            elseif p_r[3] > 0
                d_r = Vec3{_T}(zero(_T), zero(_T), H / 2)
            else
                d_r = Vec3{_T}(zero(_T), zero(_T), zero(_T))
            end
            return rotate_point(d_r, RotMat')
        end
        return Vec3{_T}(zero(_T), zero(_T), zero(_T))
    end
end

"""
    inside(p::Point{2,_T}, sill::SquareDike{2,_T})

Returns `true` if `p` is inside the 2D square-dike cross section.
"""
function inside(p::Point{2, _T}, sill::SquareDike{2, _T}) where {_T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill
    p_r = rotate_point(p - Center, RotMat)
    return abs(p_r[1]) <= W / 2 && abs(p_r[2]) <= H / 2
end

"""
    inside(p::Point{3,_T}, sill::SquareDike{3,_T})

Returns `true` if `p` is inside the 3D square-dike volume.
"""
function inside(p::Point{3, _T}, sill::SquareDike{3, _T}) where {_T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill
    p_r = rotate_point(p - Center, RotMat)
    return abs(p_r[1]) <= W / 2 && abs(p_r[2]) <= W / 2 && abs(p_r[3]) <= H / 2
end

"""
    update_abstractsill(s::SquareDike; kwargs...) -> SquareDike

Return a new square dike identical to `s` but with updated parameters.
Accepted keyword arguments: `Center`, `Angle`, `W`, `H`.
"""
function update_abstractsill(s::SquareDike; kwargs...)
    params = (
        Center = UnitValue(s.Center),
        Angle  = UnitValue(s.Angle),
        W      = UnitValue(s.W),
        H      = UnitValue(s.H),
    )
    return SquareDike(; merge(params, kwargs)...)
end
