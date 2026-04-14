using Adapt
import Base.show
import GeoParams: isdimensional

export EllipticalIntrusion

struct EllipticalIntrusion{N, _T, N1, N2, U1, U2, U3} <: AbstractSill{N, _T}
    Center::GeoUnit{Point{N, _T}, U1}
    Angle::GeoUnit{Vec{N1, _T}, U2}
    W::GeoUnit{_T, U1}
    H::GeoUnit{_T, U1}
    Lengthscale::GeoUnit{_T, U1}
    BoundingBox::Tuple
    RotMat::GeoUnit{SMatrix{N, N, _T, N2}, U3}
end
Adapt.@adapt_structure EllipticalIntrusion

isdimensional(s::EllipticalIntrusion) = isdimensional(s.W)

function EllipticalIntrusion(;
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
    return EllipticalIntrusion(Cg, convert(GeoUnit, Angle), Wg, Hg, Lengthscale, BoundingBox, convert(GeoUnit, RotMat))
end

function EllipticalIntrusion(s::EllipticalIntrusion; kwargs...)
    valid = (:Center, :Angle, :W, :H)
    all(k -> k in valid, keys(kwargs)) ||
        error("Invalid keyword for EllipticalIntrusion(s; ...). Valid keys are: $(valid)")
    base = (Center=UnitValue(s.Center), Angle=UnitValue(s.Angle), W=UnitValue(s.W), H=UnitValue(s.H))
    kw = Dict{Symbol, Any}(kwargs)
    for sym in (:W, :H)
        if haskey(kw, sym) && kw[sym] isa Number && !(kw[sym] isa typeof(oneunit(getproperty(base, sym))))
            kw[sym] = kw[sym] * oneunit(getproperty(base, sym))
        end
    end
    return EllipticalIntrusion(; merge(base, (; kw...))...)
end

function show(io::IO, s::EllipticalIntrusion)
    label = isdimensional(s) ? "dimensional units" : "nondimensional"
    println(io, "Elliptical intrusion ($label):")
    println(io, "   Center      : $((s.Center.val...,).*s.Center.unit)")
    println(io, "   Angle       : $(s.Angle.val)")
    println(io, "   Width (W)   : $(UnitValue(s.W))")
    println(io, "   Thickness(H): $(UnitValue(s.H))")
    return nothing
end

area(s::EllipticalIntrusion) = π * UnitValue(s.W) / 2 * UnitValue(s.H) / 2
volume(s::EllipticalIntrusion) = (4 / 3) * π * (UnitValue(s.W) / 2) * (UnitValue(s.W) / 2) * (UnitValue(s.H) / 2)

function hostrock_displacement(sill::EllipticalIntrusion{N, _T}, p::Point{N, _T}) where {N, _T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill
    p_r = rotate_point(p - Center, RotMat)
    AR = H / W
    if N == 2
        x = p_r[1]
        z = p_r[2]
        a = sqrt(z^2 / AR^2 + x^2)
        if a == 0
            return Vec2{_T}(zero(_T), zero(_T))
        end
        a3 = a^3
        a_inject = W / 2
        Vol_inject = (4 / 3) * π * a_inject^3 * AR
        da = ((Vol_inject + (4 / 3) * π * a3 * AR) / ((4 / 3) * π * AR))^(1 / 3) - a
        d_r = Vec2{_T}(x * (da / a), z * (da / a))
        return rotate_point(d_r, RotMat')
    else
        eq = (p_r[1]^2 + p_r[2]^2) / ((W / 2)^2) + (p_r[3]^2) / ((H / 2)^2)
        if eq == 0
            return Vec3{_T}(zero(_T), zero(_T), zero(_T))
        end
        s = eq^(1 / 6) - 1
        d_r = Vec3{_T}(p_r[1] * s, p_r[2] * s, p_r[3] * s)
        return rotate_point(d_r, RotMat')
    end
end

function inside(p::Point{2, _T}, sill::EllipticalIntrusion{2, _T}; rotate::Bool=true) where {_T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill
    p_r = p - Center
    if rotate
        p_r = rotate_point(p_r, RotMat)
    end
    eq = (p_r[1]^2) / ((W / 2)^2) + (p_r[2]^2) / ((H / 2)^2)
    return eq <= 1
end

function inside(p::Point{3, _T}, sill::EllipticalIntrusion{3, _T}; rotate::Bool=true) where {_T}
    GeoParams.@unpack_val W, H, Center, RotMat = sill
    p_r = p - Center
    if rotate
        p_r = rotate_point(p_r, RotMat)
    end
    eq = (p_r[1]^2 + p_r[2]^2) / ((W / 2)^2) + (p_r[3]^2) / ((H / 2)^2)
    return eq <= 1
end

function update_abstractsill(s::EllipticalIntrusion; kwargs...)
    params = (Center=UnitValue(s.Center), Angle=UnitValue(s.Angle), W=UnitValue(s.W), H=UnitValue(s.H))
    return EllipticalIntrusion(; merge(params, kwargs)...)
end
