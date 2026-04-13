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
end
MogiSphere(args...) = MogiSphere(convert.(GeoUnit, args)...)
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
    return MogiSphere(Center, r, ΔP, G, ν)
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
end
McTigueSphere(args...) = McTigueSphere(convert.(GeoUnit, args)...)
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
    return McTigueSphere(Center, r, ΔP, G, ν)
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
function inside(p::Point{N, _T}, sill::MogiSphere{N, _T}) where {N, _T}
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
function inside(p::Point{N, _T}, sill::McTigueSphere{N, _T}) where {N, _T}
    GeoParams.@unpack_val r, Center = sill
    Δ = p - Center
    dist_sq = zero(_T)
    for i in 1:N
        dist_sq += Δ[i]^2
    end
    return dist_sq <= r^2
end
