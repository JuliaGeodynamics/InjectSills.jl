# Centre of dilatation in an elastic half-space (Mindlin & Cheng, 1950): the point-source limit of a
# pressurized sphere with the complete interior displacement field.

export MindlinSphere

"""
    MindlinSphere{N,_T}

Pressurized spherical cavity in an elastic half-space, point-source limit (centre of
dilatation), with the complete displacement field inside the half-space (Mindlin & Cheng,
1950). At the free surface it coincides with the Mogi (1958) solution, and close to the
source it tends to the full-space field `ΔV/(4π) · r/|r|³`, so the volume displaced by the
source equals its volume change `ΔV = π ΔP r³/G`. This is the solution to use when the
displacement is needed at depth, for instance to move host rock around an inflating or
draining body: [`MogiSphere`](@ref) evaluates the *surface* formula at every point, which
is exact at `z = 0` but dilates the interior by `4(1 - ν)` times the volume change.

The free surface is the plane `z = 0`, the source sits at depth `c = -Center[N]` below it.
With `d = c - z`, `R₁` the distance to the source and `R` the distance to its image above
the surface, the displacement at `(x, y, z ≤ 0)` relative to the source is

    u = ΔV/(4π) [ (x, y, c + z)/R₁³ − (x, y, d)/R³ ]
      + ΔV/(2π) · 2(1 − ν) (x, y, d)/R³
      + ΔV/(2π) · z (3 x d/R⁵, 3 y d/R⁵, (1 − 3 d²/R²)/R³),

the closed form of the inflation point source of Okada (1992, Table 2) written in terms of
the volume change. A 2D source is the vertical plane through the axis (`y = 0`), as for
[`PennyShapedSill`](@ref) and [`MogiSphere`](@ref).

Parameters:
====
- `Center::Point{N, _T}`  - Center of the source
- `r::_T`                 - Radius of the cavity
- `ΔP::_T`                - Overpressure
- `G::_T`                 - Shear modulus
- `ν::_T`                 - Poisson's ratio

References:
===
  Mindlin, R.D., Cheng, D.H. (1950): Nuclei of strain in the semi-infinite solid.
  J. Appl. Phys. 21, 926–930. https://doi.org/10.1063/1.1699785

  Okada, Y. (1992): Internal deformation due to shear and tensile faults in a half-space.
  Bull. Seismol. Soc. Am. 82, 1018–1040.
"""
struct MindlinSphere{N, _T, U1, U2, U3} <: AbstractSill{N, _T}
    Center::GeoUnit{Point{N, _T}, U1}
    r::GeoUnit{_T, U1}
    ΔP::GeoUnit{_T, U2}
    G::GeoUnit{_T, U2}
    ν::GeoUnit{_T, U3}
    Lengthscale::GeoUnit{_T, U1}
    BoundingBox::NTuple{2, GeoUnit{Point{N, _T}, U1}}
end
Adapt.@adapt_structure MindlinSphere

isdimensional(s::MindlinSphere) = isdimensional(s.G)

"""
    MindlinSphere(; Center=Point2(0.0,-5000.0)*m, r=1500.0m, ΔP=10e6Pa, G=10e9Pa, ν=0.25*NoUnits)

Construct a half-space centre of dilatation with keyword arguments. `r` must be positive,
`ν` must satisfy `-1 < ν < 0.5` and the source must lie below the free surface `z = 0`.
"""
function MindlinSphere(;
    Center = Point2(0.0, -5000.0) * m,
    r      = 1500.0m,
    ΔP     = 10e6Pa,
    G      = 10e9Pa,
    ν      = 0.25 * NoUnits,
)
    check_positive(MindlinSphere, :r, r)
    check_poisson_ratio(MindlinSphere, ν; incompressible=false)
    Cg = convert(GeoUnit, Center)
    rg = convert(GeoUnit, r)
    Lengthscale = rg
    BoundingBox = if length(Center) == 2
        unrotated_bounding_box(Cg, rg.val, rg.val)
    else
        unrotated_bounding_box(Cg, rg.val, rg.val, rg.val)
    end
    return MindlinSphere(Cg, rg, convert(GeoUnit, ΔP), convert(GeoUnit, G), convert(GeoUnit, ν), Lengthscale, BoundingBox)
end

"""
    MindlinSphere(s::MindlinSphere; kwargs...)

Create a new `MindlinSphere` from an existing one by overriding any subset of
`Center`, `r`, `ΔP`, `G`, `ν`.
"""
function MindlinSphere(s::MindlinSphere; kwargs...)
    valid = (:Center, :r, :ΔP, :G, :ν)
    all(k -> k in valid, keys(kwargs)) ||
        error("Invalid keyword for MindlinSphere(s; ...). Valid keys are: $(valid)")
    base = (Center = UnitValue(s.Center), r = UnitValue(s.r), ΔP = UnitValue(s.ΔP), G = UnitValue(s.G), ν = UnitValue(s.ν))
    kw = Dict{Symbol,Any}(kwargs)
    for sym in (:r, :ΔP, :G)
        if haskey(kw, sym) && kw[sym] isa Number && !(kw[sym] isa typeof(oneunit(getproperty(base, sym))))
            kw[sym] = kw[sym] * oneunit(getproperty(base, sym))
        end
    end
    return MindlinSphere(; merge(base, (; kw...))...)
end

function show(io::IO, s::MindlinSphere)
    label = isdimensional(s) ? "dimensional units" : "nondimensional"
    println(io, "Mindlin half-space dilatation source ($label):")
    println(io, "   Center          : $((s.Center.val...,).*s.Center.unit)")
    println(io, "   Radius          : $(UnitValue(s.r))")
    println(io, "   Overpressure    : $(UnitValue(s.ΔP))")
    println(io, "   Shear modulus   : $(UnitValue(s.G))")
    println(io, "   Poisson's ratio : $(UnitValue(s.ν))")
    return nothing
end

# Volume and area
volume(s::MindlinSphere) = (4/3) * π * UnitValue(s.r)^3   # equivalent 3D volume of the cavity
area(s::MindlinSphere)   =         π * UnitValue(s.r)^2   # 2D cross-sectional area

"""
    d = hostrock_displacement(sill::MindlinSphere{N,_T}, p::Point{N,_T})

Displacement at `p` (with `p[N] ≤ 0`) due to a centre of dilatation of volume change
`ΔV = π ΔP r³/G` at `Center` below the free surface `z = 0`. See [`MindlinSphere`](@ref)
for the expression.
"""
function hostrock_displacement(sill::MindlinSphere{N, _T}, p::Point{N, _T}) where {N, _T}
    GeoParams.@unpack_val ν, G, r, ΔP, Center = sill

    ΔV = π * ΔP * r^3 / G
    x = p[1] - Center[1]
    y = N == 3 ? p[2] - Center[2] : zero(_T)
    z = p[N]
    c = -Center[N]                                   # depth of the source below the free surface
    d1 = c + z                                       # vertical offset from the source
    d = c - z                                        # vertical offset from the image
    ρ2 = x^2 + y^2
    R1 = max(sqrt(ρ2 + d1^2), _T(1e-8))
    R = max(sqrt(ρ2 + d^2), _T(1e-8))
    k1 = ΔV / (4 * _T(π))
    k2 = ΔV / (2 * _T(π))
    f3 = k1 / R1^3 - k1 / R^3 + k2 * 2 * (1 - ν) / R^3   # common factor of the horizontal components
    ux = x * f3 + k2 * z * 3 * x * d / R^5
    uy = y * f3 + k2 * z * 3 * y * d / R^5
    uz = k1 * d1 / R1^3 - k1 * d / R^3 + k2 * 2 * (1 - ν) * d / R^3 + k2 * z * (1 - 3 * d^2 / R^2) / R^3

    if N == 2
        return Vec2{_T}(ux, uz)
    else
        return Vec3{_T}(ux, uy, uz)
    end
end

"""
    inside(p::Point{N,_T}, sill::MindlinSphere{N,_T})

Returns `true` if `p` is inside the spherical cavity.
"""
function inside(p::Point{N, _T}, sill::MindlinSphere{N, _T}; rotate::Bool=true) where {N, _T}
    GeoParams.@unpack_val r, Center = sill
    Δ = p - Center
    dist_sq = zero(_T)
    for i in 1:N
        dist_sq += Δ[i]^2
    end
    return dist_sq <= r^2
end

"""
    update_abstractsill(s::MindlinSphere; kwargs...) -> MindlinSphere

Return a new source identical to `s` but with the specified parameters updated.
Accepted keyword arguments: `Center`, `r`, `ΔP`, `G`, `ν`.
"""
function update_abstractsill(s::MindlinSphere; kwargs...)
    params = (Center = UnitValue(s.Center), r = UnitValue(s.r), ΔP = UnitValue(s.ΔP), G = UnitValue(s.G), ν = UnitValue(s.ν))
    return MindlinSphere(; merge(params, kwargs)...)
end
