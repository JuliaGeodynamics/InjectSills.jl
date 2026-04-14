using Adapt
import Base.show
import GeoParams: isdimensional

export CylindricalDikeTopAccretion, CylindricalDikeTopAccretionFullModelAdvection

struct CylindricalDikeTopAccretion{N, _T, N1, N2, U1, U2, U3} <: AbstractSill{N, _T}
    Center::GeoUnit{Point{N, _T}, U1}
    Angle::GeoUnit{Vec{N1, _T}, U2}
    W::GeoUnit{_T, U1}
    H::GeoUnit{_T, U1}
    RotMat::GeoUnit{SMatrix{N, N, _T, N2}, U3}
end
Adapt.@adapt_structure CylindricalDikeTopAccretion

struct CylindricalDikeTopAccretionFullModelAdvection{N, _T, N1, N2, U1, U2, U3} <: AbstractSill{N, _T}
    Center::GeoUnit{Point{N, _T}, U1}
    Angle::GeoUnit{Vec{N1, _T}, U2}
    W::GeoUnit{_T, U1}
    H::GeoUnit{_T, U1}
    RotMat::GeoUnit{SMatrix{N, N, _T, N2}, U3}
end
Adapt.@adapt_structure CylindricalDikeTopAccretionFullModelAdvection

isdimensional(s::CylindricalDikeTopAccretion) = isdimensional(s.W)
isdimensional(s::CylindricalDikeTopAccretionFullModelAdvection) = isdimensional(s.W)

function CylindricalDikeTopAccretion(;
    Center = Point2(0.0, -5000.0) * m,
    Angle  = Vec1(0.0) * NoUnits,
    W      = 2000.0m,
    H      = 100.0m,
)
    @assert length(Center) == length(Angle) + 1
    RotMat = RotationMatrix(ustrip.(Angle))
    return CylindricalDikeTopAccretion(
        convert(GeoUnit, Center),
        convert(GeoUnit, Angle),
        convert(GeoUnit, W),
        convert(GeoUnit, H),
        convert(GeoUnit, RotMat),
    )
end

function CylindricalDikeTopAccretionFullModelAdvection(;
    Center = Point2(0.0, -5000.0) * m,
    Angle  = Vec1(0.0) * NoUnits,
    W      = 2000.0m,
    H      = 100.0m,
)
    @assert length(Center) == length(Angle) + 1
    RotMat = RotationMatrix(ustrip.(Angle))
    return CylindricalDikeTopAccretionFullModelAdvection(
        convert(GeoUnit, Center),
        convert(GeoUnit, Angle),
        convert(GeoUnit, W),
        convert(GeoUnit, H),
        convert(GeoUnit, RotMat),
    )
end

function CylindricalDikeTopAccretion(s::CylindricalDikeTopAccretion; kwargs...)
    valid = (:Center, :Angle, :W, :H)
    all(k -> k in valid, keys(kwargs)) ||
        error("Invalid keyword for CylindricalDikeTopAccretion(s; ...). Valid keys are: $(valid)")
    base = (Center=UnitValue(s.Center), Angle=UnitValue(s.Angle), W=UnitValue(s.W), H=UnitValue(s.H))
    kw = Dict{Symbol, Any}(kwargs)
    for sym in (:W, :H)
        if haskey(kw, sym) && kw[sym] isa Number && !(kw[sym] isa typeof(oneunit(getproperty(base, sym))))
            kw[sym] = kw[sym] * oneunit(getproperty(base, sym))
        end
    end
    return CylindricalDikeTopAccretion(; merge(base, (; kw...))...)
end

function CylindricalDikeTopAccretionFullModelAdvection(s::CylindricalDikeTopAccretionFullModelAdvection; kwargs...)
    valid = (:Center, :Angle, :W, :H)
    all(k -> k in valid, keys(kwargs)) ||
        error("Invalid keyword for CylindricalDikeTopAccretionFullModelAdvection(s; ...). Valid keys are: $(valid)")
    base = (Center=UnitValue(s.Center), Angle=UnitValue(s.Angle), W=UnitValue(s.W), H=UnitValue(s.H))
    kw = Dict{Symbol, Any}(kwargs)
    for sym in (:W, :H)
        if haskey(kw, sym) && kw[sym] isa Number && !(kw[sym] isa typeof(oneunit(getproperty(base, sym))))
            kw[sym] = kw[sym] * oneunit(getproperty(base, sym))
        end
    end
    return CylindricalDikeTopAccretionFullModelAdvection(; merge(base, (; kw...))...)
end

function show(io::IO, s::CylindricalDikeTopAccretion)
    label = isdimensional(s) ? "dimensional units" : "nondimensional"
    println(io, "Cylindrical dike top-accretion ($label):")
    println(io, "   Center      : $((s.Center.val...,).*s.Center.unit)")
    println(io, "   Angle       : $(s.Angle.val)")
    println(io, "   Width (W)   : $(UnitValue(s.W))")
    println(io, "   Thickness(H): $(UnitValue(s.H))")
    return nothing
end

function show(io::IO, s::CylindricalDikeTopAccretionFullModelAdvection)
    label = isdimensional(s) ? "dimensional units" : "nondimensional"
    println(io, "Cylindrical dike top-accretion full-model advection ($label):")
    println(io, "   Center      : $((s.Center.val...,).*s.Center.unit)")
    println(io, "   Angle       : $(s.Angle.val)")
    println(io, "   Width (W)   : $(UnitValue(s.W))")
    println(io, "   Thickness(H): $(UnitValue(s.H))")
    return nothing
end

area(s::CylindricalDikeTopAccretion) = UnitValue(s.W) * UnitValue(s.H)
volume(s::CylindricalDikeTopAccretion) = π * (UnitValue(s.W) / 2)^2 * UnitValue(s.H)
area(s::CylindricalDikeTopAccretionFullModelAdvection) = UnitValue(s.W) * UnitValue(s.H)
volume(s::CylindricalDikeTopAccretionFullModelAdvection) = π * (UnitValue(s.W) / 2)^2 * UnitValue(s.H)

function hostrock_displacement(sill::CylindricalDikeTopAccretion{N, _T}, p::Point{N, _T}) where {N, _T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill
    p_r = rotate_point(p - Center, RotMat)
    if N == 2
        if p_r[2] <= 0 && p_r[1] <= W / 2
            return rotate_point(Vec2{_T}(zero(_T), -H), RotMat')
        end
        return Vec2{_T}(zero(_T), zero(_T))
    else
        if p_r[3] <= 0 && (p_r[1]^2 + p_r[2]^2) <= (W / 2)^2
            return rotate_point(Vec3{_T}(zero(_T), zero(_T), -H), RotMat')
        end
        return Vec3{_T}(zero(_T), zero(_T), zero(_T))
    end
end

function hostrock_displacement(sill::CylindricalDikeTopAccretionFullModelAdvection{N, _T}, p::Point{N, _T}) where {N, _T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill
    p_r = rotate_point(p - Center, RotMat)
    if N == 2
        if p_r[2] <= 0
            return rotate_point(Vec2{_T}(zero(_T), -H), RotMat')
        end
        return Vec2{_T}(zero(_T), zero(_T))
    else
        if p_r[3] <= 0 && (p_r[1]^2 + p_r[2]^2) <= (W / 2)^2
            return rotate_point(Vec3{_T}(zero(_T), zero(_T), -H), RotMat')
        end
        return Vec3{_T}(zero(_T), zero(_T), zero(_T))
    end
end

function inside(p::Point{2, _T}, sill::CylindricalDikeTopAccretion{2, _T}) where {_T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill
    p_r = rotate_point(p - Center, RotMat)
    return abs(p_r[1]) <= W / 2 && abs(p_r[2]) <= H / 2
end
function inside(p::Point{3, _T}, sill::CylindricalDikeTopAccretion{3, _T}) where {_T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill
    p_r = rotate_point(p - Center, RotMat)
    return abs(p_r[1]) <= W / 2 && abs(p_r[2]) <= W / 2 && abs(p_r[3]) <= H / 2
end

function inside(p::Point{2, _T}, sill::CylindricalDikeTopAccretionFullModelAdvection{2, _T}) where {_T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill
    p_r = rotate_point(p - Center, RotMat)
    return abs(p_r[1]) <= W / 2 && abs(p_r[2]) <= H / 2
end
function inside(p::Point{3, _T}, sill::CylindricalDikeTopAccretionFullModelAdvection{3, _T}) where {_T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill
    p_r = rotate_point(p - Center, RotMat)
    return abs(p_r[1]) <= W / 2 && abs(p_r[2]) <= W / 2 && abs(p_r[3]) <= H / 2
end

function update_abstractsill(s::CylindricalDikeTopAccretion; kwargs...)
    params = (Center=UnitValue(s.Center), Angle=UnitValue(s.Angle), W=UnitValue(s.W), H=UnitValue(s.H))
    return CylindricalDikeTopAccretion(; merge(params, kwargs)...)
end

function update_abstractsill(s::CylindricalDikeTopAccretionFullModelAdvection; kwargs...)
    params = (Center=UnitValue(s.Center), Angle=UnitValue(s.Angle), W=UnitValue(s.W), H=UnitValue(s.H))
    return CylindricalDikeTopAccretionFullModelAdvection(; merge(params, kwargs)...)
end
