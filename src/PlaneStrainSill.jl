# Pressurized crack in plane strain: a sill in a Cartesian 2D model
using Adapt
import Base.show
import GeoParams: isdimensional

export PlaneStrainSill

"""
    PlaneStrainSill{_T}

Pressurized crack in plane strain in a homogeneous elastic full space (Westergaard, 1939;
Pollard & Segall, 1987): a sill in a Cartesian 2D model, i.e. of infinite extent normal to
the model plane. It is 2D only; for a 3D sill or an axisymmetric 2D model, use
[`PennyShapedSill`](@ref).

Parameters:
====
- `Center::Point{2, _T}`  - Center of sill
- `Angle::Vec{1, _T}`     - Dip angle of sill w.r.t. horizontal
- `E::_T`                 - Young's modulus
- `ν::_T`                 - Poisson's ratio
- `ΔP::_T`                - Overpressure within sill
- `Q::_T = π*R*H/2`       - Cross-sectional area of sill [m^2]
- `R::_T`                 - Half-length of sill
- `H::_T = 4*(1-ν^2)*ΔP*R/E` - Maximum opening of sill

The opening at distance `x` from the center is `H√(1 - x²/R²)`, and half of the area `Q`
crosses every line parallel to the sill on either side. `area` is `Q`, the cross-sectional
area per unit length normal to the model plane. As for the other 2D sills, `volume` is that of
the 3D body with this cross-section, the spheroid `4π R² (H/2)/3` of the penny-shaped sill.

Reference:
===
   Westergaard, H.M., 1939. Bearing pressures and cracks. J. Appl. Mech. 6, A49–A53.

   Pollard, D.D., Segall, P., 1987. Theoretical displacements and stresses near fractures in
        rock: with applications to faults, joints, veins, dikes, and solution surfaces.
        In: Atkinson, B.K. (Ed.), Fracture Mechanics of Rock. Academic Press, 277–349.
"""
struct PlaneStrainSill{_T, U1, U2, U3, U4, U5} <: AbstractSill{2, _T}
    Center::GeoUnit{Point{2, _T}, U1}
    Angle::GeoUnit{Vec{1, _T}, U2}
    E::GeoUnit{_T, U3}
    ν::GeoUnit{_T, U4}
    ΔP::GeoUnit{_T, U3}
    Q::GeoUnit{_T, U5}
    R::GeoUnit{_T, U1}
    H::GeoUnit{_T, U1}
    Lengthscale::GeoUnit{_T, U1}
    BoundingBox::NTuple{2, GeoUnit{Point{2, _T}, U1}}
    RotMat::GeoUnit{SMatrix{2, 2, _T, 4}, U4}
end
Adapt.@adapt_structure PlaneStrainSill

isdimensional(s::PlaneStrainSill) = isdimensional(s.E)

"""
    PlaneStrainSill(; R=nothing, Q=nothing, ΔP=nothing, H=nothing, E=1.5e10Pa, ν=0.3*NoUnits, Angle=Vec1(0.0)*NoUnits, Center=Point2(0.0)*m)

Defines a pressurized crack in plane strain. `Center` must be a 2D point. The geometry is
set by at most two of `R`, `H`, `ΔP`, `Q`, related by `H = 4(1-ν²) ΔP R/E` and
`Q = π R H/2`:
- any two of them (`R`+`H`, `R`+`Q`, `R`+`ΔP`, `H`+`Q`, `H`+`ΔP`, `Q`+`ΔP`)
- only `Q` (with `ΔP = 1e6Pa`) or only `ΔP` (with the `R` of the default sill)
- none of them (defaults `ΔP = 1e6Pa`, `Q = 1000m^2`)

Giving three or four of them, or only `R` or only `H`, throws an `ArgumentError`.
`R`, `H`, `ΔP` and `Q` must be positive, and `ν` must satisfy `-1 < ν ≤ 0.5`.
"""
function PlaneStrainSill(; R=nothing, Q=nothing, ΔP=nothing, H=nothing, E=1.5e10Pa, ν=0.3*NoUnits, Angle=Vec1(0.0)*NoUnits, Center=Point2(0.0)*m)
    length(Center) == 2 ||
        throw(ArgumentError("PlaneStrainSill: `Center` must be 2D; got $(length(Center)) coordinates. A plane-strain sill exists in 2D only; use PennyShapedSill in 3D"))
    length(Angle) == 1 ||
        throw(ArgumentError("PlaneStrainSill: `Angle` must have one component (the dip); got $(length(Angle))"))
    check_poisson_ratio(PlaneStrainSill, ν; incompressible=true)
    check_geometry(PlaneStrainSill, R, H, ΔP, Q)

    k = 4 * (1 - ν^2) / E     # H = k ΔP R
    if isnothing(R) && isnothing(H)
        if isnothing(Q) && !isnothing(ΔP)
            # the half-length of the default sill
            R = sqrt(2 * 1000.0m^2 / (k * 1e6Pa * π))
            H = k * ΔP * R
            Q = R * H * π / 2
        else
            ΔP = isnothing(ΔP) ? 1e6Pa : ΔP
            Q  = isnothing(Q) ? 1000.0m^2 : Q
            R  = sqrt(2 * Q / (k * ΔP * π))
            H  = k * ΔP * R
        end
    elseif !isnothing(R) && !isnothing(H)
        ΔP = H / (k * R)
        Q  = R * H * π / 2
    elseif !isnothing(R) && !isnothing(Q)
        H  = 2 * Q / (R * π)
        ΔP = H / (k * R)
    elseif !isnothing(R) && !isnothing(ΔP)
        H  = k * ΔP * R
        Q  = R * H * π / 2
    elseif !isnothing(H) && !isnothing(Q)
        R  = 2 * Q / (H * π)
        ΔP = H / (k * R)
    elseif !isnothing(H) && !isnothing(ΔP)
        R  = H / (k * ΔP)
        Q  = R * H * π / 2
    end

    RotMat = RotationMatrix(ustrip.(Angle))
    Cg = convert(GeoUnit, Center)
    Rg = convert(GeoUnit, R)
    Hg = convert(GeoUnit, H)
    BoundingBox = unrotated_bounding_box(Cg, Rg.val, Hg.val / 2)
    return PlaneStrainSill(Cg, convert(GeoUnit, Angle), convert(GeoUnit, E), convert(GeoUnit, ν),
                           convert(GeoUnit, ΔP), convert(GeoUnit, Q), Rg, Hg, Hg, BoundingBox,
                           convert(GeoUnit, RotMat))
end

"""
    PlaneStrainSill(s::PlaneStrainSill; kwargs...)

Create a new plane-strain sill from an existing one by changing any of `Center`, `Angle`,
`E`, `ν`, `R`, `H`, `ΔP`, `Q`, as [`update_abstractsill`](@ref) does. Numbers without unit
for `E`, `ΔP`, `Q`, `R` or `H` take the unit of the value they replace.
"""
PlaneStrainSill(s::PlaneStrainSill; kwargs...) = copy_RHΔPQ(PlaneStrainSill, s; kwargs...)

function show(io::IO, g::PlaneStrainSill)
    label = isdimensional(g) ? "dimensional units" : "nondimensional"
    println(io, "Plane-strain sill ($label):")
    println(io, "   Center                  : $((g.Center.val...,).*g.Center.unit) ")
    println(io, "   Angle [degree]          : $((g.Angle.val...,)) ")
    println(io, "   Young's modulus         : $(UnitValue(g.E)) ")
    println(io, "   Poisson's ratio         : $(UnitValue(g.ν)) ")
    println(io, "   Overpressure            : $(UnitValue(g.ΔP)) ")
    println(io, "   Cross-sectional area    : $(UnitValue(g.Q)) ")
    println(io, "   Maximum opening         : $(UnitValue(g.H)) ")
    println(io, "   Half-length             : $(UnitValue(g.R)) ")
    return nothing
end

area(s::PlaneStrainSill) = UnitValue(s.Q)
volume(s::PlaneStrainSill) = 4 / 3 * π * UnitValue(s.R)^2 * (UnitValue(s.H) / 2)

"""
    d = hostrock_displacement(sill::PlaneStrainSill{_T}, p::Point{2, _T})

Host rock displacement at `p` caused by the opening of a plane-strain sill.
"""
function hostrock_displacement(sill::PlaneStrainSill{_T}, p::Point{2, _T}) where {_T}
    # Sill frame: x along the sill, y normal to it, half-length a = R, internal pressure p.
    # With ζ = x + i y and s = √(ζ-a)√(ζ+a) (branch cut on the crack, s ~ ζ at infinity),
    # Ẑ = -p a²/(ζ + s) and Z = dẐ/dζ = p a²/(s(ζ + s)),
    #     2G u_x = (1-2ν) Re Ẑ - y Im Z,    2G u_y = 2(1-ν) Im Ẑ - y Re Z    (y ≥ 0),
    # and u_x is even, u_y odd in y. The opening u_y(x,0⁺) - u_y(x,0⁻) = 4(1-ν) p a/(2G) √(1-x²/a²)
    # is H√(1-x²/R²) for p/(2G) = H/(4(1-ν)R). The kernel works with ξ = ζ/R and
    # w = 1/(ξ + s/R), so that Ẑ = -pR w and Z = p w/(s/R) are formed without dimensional powers.
    GeoParams.@unpack_val ν, R, H, Center, RotMat = sill
    Δ  = rotate_point(p - Center, RotMat)
    y  = abs(Δ[2]) / R
    ξ  = complex(Δ[1] / R, y)
    s  = sqrt(ξ - 1) * sqrt(ξ + 1)
    w  = cdiv(one(_T), ξ + s)
    # y·Z vanishes on the sill plane; Z itself is infinite at the tips (s = 0).
    yZ = y > 0 ? y * cdiv(w, s) : zero(w)
    k  = H / (4 * (1 - ν))
    ux = -k * ((1 - 2ν) * real(w) + imag(yZ))
    uy = -k * (2 * (1 - ν) * imag(w) + real(yZ))
    d  = Vec2{_T}(ux, Δ[2] < 0 ? -uy : uy)
    return rotate_point(d, RotMat')
end

"""
    inside(p::Point{2, _T}, sill::PlaneStrainSill{_T})

Checks if the point `p` is inside the sill.
"""
function inside(p::Point{2, _T}, sill::PlaneStrainSill{_T}; rotate::Bool=true) where {_T}
    GeoParams.@unpack_val R, H, Center, RotMat = sill
    p_r = p - Center
    if rotate
        p_r = rotate_point(p_r, RotMat)
    end
    return (p_r[1] / R)^2 + (p_r[2] / (H / 2))^2 <= 1
end

"""
    update_abstractsill(s::PlaneStrainSill; kwargs...) -> PlaneStrainSill

Return a new sill identical to `s` but with the specified parameters updated.
Accepted keyword arguments: `Center`, `Angle`, `E`, `ν`, `R`, `H`, `ΔP`, `Q`.

When a single geometric parameter is changed, the remaining ones are kept
consistent using these defaults:
- only `R` → keep `H`, recompute `ΔP` and `Q`
- only `H` → keep `R`, recompute `ΔP` and `Q`
- only `ΔP` → keep `H`, recompute `R` and `Q`
- only `Q` → keep `H`, recompute `R` and `ΔP`
"""
update_abstractsill(s::PlaneStrainSill; kwargs...) = update_RHΔPQ(PlaneStrainSill, s; kwargs...)
