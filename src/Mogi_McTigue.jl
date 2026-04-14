using Adapt
import Base.show
import GeoParams: isdimensional

export MogiSphere, McTigueSphere

# ---- MogiSphere ----------------------------------------------------------

"""
    MogiSphere{N,_T}

Pressurized spherical cavity in an elastic half-space (Mogi, 1958).

Parameters:
===
- `Center::Point{N,_T}` - centre of sphere
- `r::_T`               - radius [m]
- `ΔP::_T`              - overpressure [Pa]
- `G::_T`               - shear modulus [Pa]
- `ν::_T`               - Poisson's ratio [-]

Reference:
===
  Mogi, K. (1958): Relations between the eruptions of various volcanoes and the
  deformations of the ground surfaces around them.
  Bull. Earthq. Res. Inst. Univ. Tokyo, 36, 99–134.
"""
struct MogiSphere{N, _T, U1, U2, U3} <: AbstractSill{N, _T}
    Center::GeoUnit{Point{N, _T}, U1}
    r::GeoUnit{_T, U1}
    ΔP::GeoUnit{_T, U2}
    G::GeoUnit{_T, U2}
    ν::GeoUnit{_T, U3}
    Lengthscale::GeoUnit{_T, U1}
    BoundingBox::Tuple
end
Adapt.@adapt_structure MogiSphere

isdimensional(s::MogiSphere) = isdimensional(s.G)

"""
    MogiSphere(; Center=Point2(0.0,-5000.0)*m, r=1500.0m, ΔP=10e6Pa, G=10e9Pa, ν=0.25*NoUnits)

Construct a Mogi spherical pressure source with keyword arguments.
"""
function MogiSphere(;
    Center = Point2(0.0, -5000.0) * m,
    r      = 1500.0m,
    ΔP     = 10e6Pa,
    G      = 10e9Pa,
    ν      = 0.25 * NoUnits,
)
    Cg = convert(GeoUnit, Center)
    rg = convert(GeoUnit, r)
    Lengthscale = rg
    BoundingBox = if length(Center) == 2
        unrotated_bounding_box(Cg, rg.val, rg.val)
    else
        unrotated_bounding_box(Cg, rg.val, rg.val, rg.val)
    end
    return MogiSphere(Cg, rg, convert(GeoUnit, ΔP), convert(GeoUnit, G), convert(GeoUnit, ν), Lengthscale, BoundingBox)
end

"""
    MogiSphere(s::MogiSphere; kwargs...)

Create a new `MogiSphere` from an existing one by overriding any subset of
`Center`, `r`, `ΔP`, `G`, `ν`.
"""
function MogiSphere(s::MogiSphere; kwargs...)
    valid = (:Center, :r, :ΔP, :G, :ν)
    all(k -> k in valid, keys(kwargs)) ||
        error("Invalid keyword for MogiSphere(s; ...). Valid keys are: $(valid)")

    base = (
        Center = UnitValue(s.Center),
        r      = UnitValue(s.r),
        ΔP     = UnitValue(s.ΔP),
        G      = UnitValue(s.G),
        ν      = UnitValue(s.ν),
    )

    kw = Dict{Symbol,Any}(kwargs)
    for sym in (:r, :ΔP, :G)
        if haskey(kw, sym) && kw[sym] isa Number && !(kw[sym] isa typeof(oneunit(getproperty(base, sym))))
            kw[sym] = kw[sym] * oneunit(getproperty(base, sym))
        end
    end

    return MogiSphere(; merge(base, (; kw...))...)
end

function show(io::IO, s::MogiSphere)
    label = isdimensional(s) ? "dimensional units" : "nondimensional"
    println(io, "Mogi sphere ($label):")
    println(io, "   Center          : $((s.Center.val...,).*s.Center.unit)")
    println(io, "   Radius          : $(UnitValue(s.r))")
    println(io, "   Overpressure    : $(UnitValue(s.ΔP))")
    println(io, "   Shear modulus   : $(UnitValue(s.G))")
    println(io, "   Poisson's ratio : $(UnitValue(s.ν))")
    return nothing
end

# Volume and area
volume(s::MogiSphere)   =  (4/3) * π * UnitValue(s.r)^3   # equivalent 3D volume
area(s::MogiSphere)     =          π * UnitValue(s.r)^2   # 2D cross-sectional area


# ---- McTigueSphere -------------------------------------------------------

"""
    McTigueSphere{N,_T}

Pressurized spherical cavity in an elastic half-space, higher-order solution
(McTigue, 1987). Same parameters as `MogiSphere`; adds a finite-radius
correction proportional to `(r/d)³`.

Reference:
===
  McTigue, D.F. (1987): Elastic stress and deformation near a finite spherical
  magma body: resolution of the point source paradox.
  J. Geophys. Res. 92 (B12), 12931–12940.
"""
struct McTigueSphere{N, _T, U1, U2, U3} <: AbstractSill{N, _T}
    Center::GeoUnit{Point{N, _T}, U1}
    r::GeoUnit{_T, U1}
    ΔP::GeoUnit{_T, U2}
    G::GeoUnit{_T, U2}
    ν::GeoUnit{_T, U3}
    Lengthscale::GeoUnit{_T, U1}
    BoundingBox::Tuple
end
Adapt.@adapt_structure McTigueSphere

isdimensional(s::McTigueSphere) = isdimensional(s.G)

"""
    McTigueSphere(; Center=Point2(0.0,-5000.0)*m, r=1500.0m, ΔP=10e6Pa, G=10e9Pa, ν=0.25*NoUnits)
"""
function McTigueSphere(;
    Center = Point2(0.0, -5000.0) * m,
    r      = 1500.0m,
    ΔP     = 10e6Pa,
    G      = 10e9Pa,
    ν      = 0.25 * NoUnits,
)
    Cg = convert(GeoUnit, Center)
    rg = convert(GeoUnit, r)
    Lengthscale = rg
    BoundingBox = if length(Center) == 2
        unrotated_bounding_box(Cg, rg.val, rg.val)
    else
        unrotated_bounding_box(Cg, rg.val, rg.val, rg.val)
    end
    return McTigueSphere(Cg, rg, convert(GeoUnit, ΔP), convert(GeoUnit, G), convert(GeoUnit, ν), Lengthscale, BoundingBox)
end

"""
    McTigueSphere(s::McTigueSphere; kwargs...)

Create a new `McTigueSphere` from an existing one by overriding any subset of
`Center`, `r`, `ΔP`, `G`, `ν`.
"""
function McTigueSphere(s::McTigueSphere; kwargs...)
    valid = (:Center, :r, :ΔP, :G, :ν)
    all(k -> k in valid, keys(kwargs)) ||
        error("Invalid keyword for McTigueSphere(s; ...). Valid keys are: $(valid)")

    base = (
        Center = UnitValue(s.Center),
        r      = UnitValue(s.r),
        ΔP     = UnitValue(s.ΔP),
        G      = UnitValue(s.G),
        ν      = UnitValue(s.ν),
    )

    kw = Dict{Symbol,Any}(kwargs)
    for sym in (:r, :ΔP, :G)
        if haskey(kw, sym) && kw[sym] isa Number && !(kw[sym] isa typeof(oneunit(getproperty(base, sym))))
            kw[sym] = kw[sym] * oneunit(getproperty(base, sym))
        end
    end

    return McTigueSphere(; merge(base, (; kw...))...)
end

function show(io::IO, s::McTigueSphere)
    label = isdimensional(s) ? "dimensional units" : "nondimensional"
    println(io, "McTigue sphere ($label):")
    println(io, "   Center          : $((s.Center.val...,).*s.Center.unit)")
    println(io, "   Radius          : $(UnitValue(s.r))")
    println(io, "   Overpressure    : $(UnitValue(s.ΔP))")
    println(io, "   Shear modulus   : $(UnitValue(s.G))")
    println(io, "   Poisson's ratio : $(UnitValue(s.ν))")
    return nothing
end

# Volume and area
volume(s::McTigueSphere)   =  (4/3) * π * UnitValue(s.r)^3   # equivalent 3D volume
area(s::McTigueSphere)     =          π * UnitValue(s.r)^2   # 2D cross-sectional area

# ---- hostrock_displacement -----------------------------------------------

"""
    d = hostrock_displacement(sill::MogiSphere{N,_T}, p::Point{N,_T})

Displacement at `p` due to a Mogi spherical pressure source.
The displacement vector is simply proportional to `(p - Center)`:

    U = (r³ ΔP (1−ν) / G) · (p − Center) / |p − Center|³
"""
function hostrock_displacement(sill::MogiSphere{N, _T}, p::Point{N, _T}) where {N, _T}
    GeoParams.@unpack_val ν, G, r, ΔP, Center = sill

    Δ = p - Center
    R_sq = zero(_T)
    for i in 1:N
        R_sq += Δ[i]^2
    end
    R = sqrt(R_sq)
    if R < 1e-8; R = convert(_T, 1e-8); end

    C = r^3 * ΔP * (1 - ν) / (G * R^3)

    if N == 2
        return Vec2{_T}(C * Δ[1], C * Δ[2])
    else
        return Vec3{_T}(C * Δ[1], C * Δ[2], C * Δ[3])
    end
end

"""
    d = hostrock_displacement(sill::McTigueSphere{N,_T}, p::Point{N,_T})

Displacement at `p` due to a McTigue spherical pressure source.
Adds a `(r/d)³` finite-size correction to the Mogi solution, where `d`
is the vertical distance from source centre to observation point.
"""
function hostrock_displacement(sill::McTigueSphere{N, _T}, p::Point{N, _T}) where {N, _T}
    GeoParams.@unpack_val ν, G, r, ΔP, Center = sill

    Δ = p - Center
    R_sq = zero(_T)
    for i in 1:N
        R_sq += Δ[i]^2
    end
    R = sqrt(R_sq)
    if R < 1e-8; R = convert(_T, 1e-8); end

    d = abs(Δ[N])
    if d < 1e-8; d = convert(_T, 1e-8); end

    C = r^3 * ΔP * (1 - ν) / (G * R^3)

    # McTigue finite-source correction
    term1 = (r / d)^3
    term2 = (1 + ν) / (2 * (-7 + 5ν))
    term3 = 15 * d^2 * (-2 + ν) / (4 * R_sq * (-7 + 5ν))
    C    *= 1 + term1 * (term2 + term3)

    if N == 2
        return Vec2{_T}(C * Δ[1], C * Δ[2])
    else
        return Vec3{_T}(C * Δ[1], C * Δ[2], C * Δ[3])
    end
end


# ---- inside --------------------------------------------------------------

"""
    inside(p::Point{N,_T}, sill::MogiSphere{N,_T})

Returns `true` if `p` is inside the spherical cavity.
"""
function inside(p::Point{N, _T}, sill::MogiSphere{N, _T}; rotate::Bool=true) where {N, _T}
    GeoParams.@unpack_val r, Center = sill
    Δ = p - Center
    dist_sq = zero(_T)
    for i in 1:N
        dist_sq += Δ[i]^2
    end
    return dist_sq <= r^2
end

"""
    inside(p::Point{N,_T}, sill::McTigueSphere{N,_T})

Returns `true` if `p` is inside the spherical cavity.
"""
function inside(p::Point{N, _T}, sill::McTigueSphere{N, _T}; rotate::Bool=true) where {N, _T}
    GeoParams.@unpack_val r, Center = sill
    Δ = p - Center
    dist_sq = zero(_T)
    for i in 1:N
        dist_sq += Δ[i]^2
    end
    return dist_sq <= r^2
end


"""
    update_abstractsill(s::MogiSphere; kwargs...) -> MogiSphere
    update_abstractsill(s::McTigueSphere; kwargs...) -> McTigueSphere

Return a new sphere identical to `s` but with the specified parameters updated.
Accepted keyword arguments: `Center`, `r`, `ΔP`, `G`, `ν`.

# Example
```julia
s2 = update_abstractsill(s, r = 2000.0m, ΔP = 15e6Pa)
```
"""
function update_abstractsill(s::MogiSphere; kwargs...)
    params = (
        Center = UnitValue(s.Center),
        r      = UnitValue(s.r),
        ΔP     = UnitValue(s.ΔP),
        G      = UnitValue(s.G),
        ν      = UnitValue(s.ν),
    )
    return MogiSphere(; merge(params, kwargs)...)
end

function update_abstractsill(s::McTigueSphere; kwargs...)
    params = (
        Center = UnitValue(s.Center),
        r      = UnitValue(s.r),
        ΔP     = UnitValue(s.ΔP),
        G      = UnitValue(s.G),
        ν      = UnitValue(s.ν),
    )
    return McTigueSphere(; merge(params, kwargs)...)
end
