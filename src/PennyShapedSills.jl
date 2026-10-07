# Penny shapes sills embedded in a homogeneous elastic halfspace
using Adapt
import Base.show
import GeoParams: isdimensional

export PennyShapedSill, set_penny_shaped_sill, hostrock_displacement
export inside
export update_abstractsill


"""
    PennyShapedSill{N,_T}

Holds information about a penny shaped sill in 2D or 3D.

Parameters:
====
- `Center::Point{N, _T}`  - Center of sill
- `Angle::Vec{N1, _T}`    - Dip and strike (in 3D) angle of sill w.r.t. horizontal
- `E::_T`                 - Young's modulus
- `ν::_T`                 - Poisson's ratio
- `ΔP::_T`                - Overpressure within sill
- `Q::_T`                 - Total injected volume of sill [m^3]
- `R::_T = cbrt(3*E*Q/(16*(1-ν^2)*ΔP))`     - Radius (3D) or half-length of the axisymmetric section (2D)
- `H::_T = 8*(1-ν^2)*ΔP*R/(π*E)`                  - Maximum thickness of sill

A 2D sill is the axisymmetric section of the penny: its displacement is the 3D solution in
the vertical plane through the sill axis, and `ΔP` and `Q` are those of the 3D penny. For a
sill in a Cartesian 2D (plane-strain) model, use [`PlaneStrainSill`](@ref).

Reference:
===
   Sun, R.J., 1969. Theoretical size of hydraulically induced horizontal fractures and
        corresponding surface uplift in an idealized medium. J. Geophys. Res. 74, 5995–6011.
        https://doi.org/10.1029/JB074i025p05995
"""
struct PennyShapedSill{N, _T, N1, N2, U1, U2, U3, U4, U5} <: AbstractSill{N,_T}
    Center::GeoUnit{Point{N, _T},U1}  # m
    Angle::GeoUnit{Vec{N1, _T},U2}    # degrees
    E::GeoUnit{_T,U3}   # in Pa
    ν::GeoUnit{_T,U4}   # []
    ΔP::GeoUnit{_T,U3}  # Pa
    Q::GeoUnit{_T,U5}   # m^3
    R::GeoUnit{_T,U1}   # m
    H::GeoUnit{_T,U1}   # m
    Lengthscale::GeoUnit{_T,U1}
    BoundingBox::NTuple{2, GeoUnit{Point{N, _T}, U1}}
    RotMat::GeoUnit{SMatrix{N,N,_T,N2},U4}             # rotation matrix (precomputed for efficiency)
    RotMat_negative::GeoUnit{SMatrix{N,N,_T,N2},U4}    # with negative angle
end
Adapt.@adapt_structure PennyShapedSill

isdimensional(s::PennyShapedSill) = isdimensional(s.E)

"""
    PennyShapedSill(; R=nothing,  Q=nothing, ΔP=nothing, H=nothing, E=1.5e10Pa, ν=0.3*NoUnits, Angle=Vec1(0.0)*Pas, Center=Point2(0.0)*m)

Defines parameters for a penny shaped sill in a homogeneous elastic medium.
The geometry is set by at most two of `R`, `H`, `ΔP`, `Q`:
- any two of them (`R`+`H`, `R`+`Q`, `R`+`ΔP`, `H`+`Q`, `H`+`ΔP`, `Q`+`ΔP`)
- only `Q` (with `ΔP = 1e6Pa`) or only `ΔP` (with the `R` of the default sill)
- none of them (defaults `ΔP = 1e6Pa`, `Q = 1000m^3`)

Giving three or four of them, or only `R` or only `H`, throws an `ArgumentError`.
`R`, `H`, `ΔP` and `Q` must be positive, and `ν` must satisfy `-1 < ν < 0.5`
(the solution of Sun (1969) contains `1/(1-2ν)`).
"""
function PennyShapedSill(; R=nothing,  Q=nothing, ΔP=nothing, H=nothing, E=1.5e10Pa, ν=0.3*NoUnits, Angle=Vec1(0.0)*Pas, Center=Point2(0.0)*m, W=nothing)
    reject_W(PennyShapedSill, "radius", W)
    @assert length(Center)==length(Angle)+1
    check_poisson_ratio(PennyShapedSill, ν; incompressible=false)
    given = [name for (name, x) in ((:R, R), (:H, H), (:ΔP, ΔP), (:Q, Q)) if !isnothing(x)]
    if length(given) > 2 || given == [:R] || given == [:H]
        throw(ArgumentError("PennyShapedSill: got $(join(given, ", ")); specify at most two of R, H, ΔP, Q (any pair), only Q, only ΔP, or none of them"))
    end
    for (name, x) in ((:R, R), (:H, H), (:ΔP, ΔP), (:Q, Q))
        isnothing(x) || check_positive(PennyShapedSill, name, x)
    end

    if isnothing(R) && isnothing(Q) && isnothing(ΔP) && isnothing(H)
        ΔP  =  1e6*Pa
        Q   =  1000.0*m^3
        R   =  cbrt(3*E*Q/(16*(1-ν^2)*ΔP))
        H   =  8*(1-ν^2)*ΔP*R/(π*E)

    elseif isnothing(R) && isnothing(Q) && !isnothing(ΔP) && isnothing(H)
        # Use the default reference radius and recompute H and Q for the
        # user-provided overpressure.
        Q_ref = 1000.0*m^3
        ΔP_ref = 1e6*Pa
        R   =  cbrt(3*E*Q_ref/(16*(1-ν^2)*ΔP_ref))
        H   =  8*(1-ν^2)*ΔP*R/(π*E)
        Q   =  16*(1-ν^2)*ΔP*R^3/(3*E)

    elseif isnothing(R) && !isnothing(Q) && !isnothing(ΔP) && isnothing(H)
        R   =  cbrt(3*E*Q/(16*(1-ν^2)*ΔP))
        H   =  8*(1-ν^2)*ΔP*R/(π*E)

    elseif isnothing(R) && !isnothing(Q) && isnothing(ΔP) && isnothing(H)
        ΔP  =  1e6*Pa
        R   =  cbrt(3*E*Q/(16*(1-ν^2)*ΔP))
        H   =  8*(1-ν^2)*ΔP*R/(π*E)

    elseif !isnothing(R) && !isnothing(Q)
        ΔP = 3*E*Q/(16*(1-ν^2)*R^3)
        H  = (3 * Q) / (2 * R^2 * π)

    elseif !isnothing(R) && !isnothing(H)
        ΔP = (π * E * H) / (8 * (1 - ν^2) * R)
        Q = (2 * H * R^2 * π) / 3

    elseif !isnothing(R) && !isnothing(ΔP)
        H = 8 * (1 - ν^2) * ΔP * R / (π * E)
        Q = 16 * (1 - ν^2) * ΔP * R^3 / (3 * E)

    elseif !isnothing(H) && !isnothing(Q)
        R = sqrt(3 * Q / (2 * H * π))
        ΔP = (π * E * H) / (8 * (1 - ν^2) * R)

    elseif !isnothing(H) && !isnothing(ΔP)
        R = (π * E * H) / (8 * (1 - ν^2) * ΔP)
        Q = 16 * (1 - ν^2) * ΔP * R^3 / (3 * E)
    end

    # Compute rotation matrix - as this is a relatively expensive operation, we precompute & store it in the struct
    RotMat = RotationMatrix(ustrip.(Angle))
    RotMat_negative = RotMat'

    Cg = convert(GeoUnit, Center)
    Rg = convert(GeoUnit, R)
    Hg = convert(GeoUnit, H)
    Lengthscale = Hg
    BoundingBox = if length(Center) == 2
        unrotated_bounding_box(Cg, Rg.val, Hg.val / 2)
    else
        unrotated_bounding_box(Cg, Rg.val, Rg.val, Hg.val / 2)
    end

    return PennyShapedSill(
        convert(GeoUnit, Center),
        convert(GeoUnit, Angle),
        convert(GeoUnit, E),
        convert(GeoUnit, ν),
        convert(GeoUnit, ΔP),
        convert(GeoUnit, Q),
        Rg,
        Hg,
        Lengthscale,
        BoundingBox,
        convert(GeoUnit, RotMat),
        convert(GeoUnit, RotMat_negative),
    )
end


"""
    PennyShapedSill(s::PennyShapedSill; kwargs...)

Create a new penny-shaped sill from an existing one by changing any number of
parameters.

For updates of `R`, `H`, `ΔP`, or `Q`, a valid constructor combination is
selected automatically.
"""
function PennyShapedSill(s::PennyShapedSill; kwargs...)
    reject_W(PennyShapedSill, "radius", get(kwargs, :W, nothing))
    valid = (:Center, :Angle, :E, :ν, :R, :H, :ΔP, :Q)
    all(k -> k in valid, keys(kwargs)) ||
        error("Invalid keyword for PennyShapedSill(s; ...). Valid keys are: $(valid)")

    base = (
        Center = UnitValue(s.Center),
        Angle  = UnitValue(s.Angle),
        E      = UnitValue(s.E),
        ν      = UnitValue(s.ν),
        R      = UnitValue(s.R),
        H      = UnitValue(s.H),
        ΔP     = UnitValue(s.ΔP),
        Q      = UnitValue(s.Q),
    )

    kw = Dict{Symbol,Any}(kwargs)
    for sym in (:E, :ΔP, :Q, :R, :H)
        if haskey(kw, sym) && kw[sym] isa Number && !(kw[sym] isa typeof(oneunit(getproperty(base, sym))))
            kw[sym] = kw[sym] * oneunit(getproperty(base, sym))
        end
    end

    updated = merge(base, (; kw...))
    gkeys = Set{Symbol}(k for k in keys(kwargs) if k in (:R, :H, :ΔP, :Q))

    # Select the correct constructor branch for physically consistent updates.
    R  = nothing
    H  = nothing
    ΔP = nothing
    Q  = nothing

    if isempty(gkeys)
        R = updated.R
        H = updated.H
    elseif gkeys == Set([:R])
        R = updated.R
        H = updated.H
    elseif gkeys == Set([:H])
        R = updated.R
        H = updated.H
    elseif gkeys == Set([:ΔP])
        H = updated.H
        ΔP = updated.ΔP
    elseif gkeys == Set([:Q])
        H = updated.H
        Q = updated.Q
    elseif gkeys == Set([:R, :H])
        R = updated.R
        H = updated.H
    elseif gkeys == Set([:R, :Q])
        R = updated.R
        Q = updated.Q
    elseif gkeys == Set([:R, :ΔP])
        R = updated.R
        ΔP = updated.ΔP
    elseif gkeys == Set([:H, :Q])
        H = updated.H
        Q = updated.Q
    elseif gkeys == Set([:H, :ΔP])
        H = updated.H
        ΔP = updated.ΔP
    else
        error("Invalid R/H/ΔP/Q combination. Use one of: R, H, ΔP, Q, R+H, R+Q, R+ΔP, H+Q, H+ΔP")
    end

    return PennyShapedSill(; Center=updated.Center, Angle=updated.Angle, E=updated.E, ν=updated.ν, R, H, ΔP, Q)
end

# Print info in the REPL
function show(io::IO, g::PennyShapedSill)

    if isdimensional(g)
        println(io, "Penny-shaped sill in dimensional units:")
        println(io, "   Center                  : $((g.Center.val...,).*g.Center.unit) ")
    else
        println(io, "Penny-shaped sill in nondimensional units:")
        println(io, "   Center                  : $((g.Center.val...,).*g.Center.unit) ")
    end
    println(io, "   Angle [degree]          : $((g.Angle.val...,)) ")
    println(io, "   Young's modulus         : $(UnitValue(g.E)) ")
    println(io, "   Poison's ratio          : $(UnitValue(g.ν)) ")
    println(io, "   Overpressure            : $(UnitValue(g.ΔP)) ")
    println(io, "   Sill volume             : $(UnitValue(g.Q)) ")
    println(io, "   Maximum sill thickness  : $(UnitValue(g.H)) ")
    println(io, "   Radius                  : $(UnitValue(g.R)) ")

    return nothing
end

# Volume and area.
volume(s::PennyShapedSill) =  4/3*π*UnitValue(s.R)*UnitValue(s.R)*(UnitValue(s.H)/2)   #   (equivalent 3D volume = injected Q, in m^3)
area(s::PennyShapedSill) = π*UnitValue(s.R)*(UnitValue(s.H)/2)                          #   (vertical cross-section, in m^2)

"""
    d = hostrock_displacement(sill::PennyShapedSill{N,_T}, p::Point{N, _T})

Host rock displacement caused by opening of a penny shaped sill at point `p`
"""
function hostrock_displacement(sill::PennyShapedSill{N,_T}, p::Point{N, _T}) where {N,_T}

    GeoParams.@unpack_val ν,E,R,H, ΔP, Center, Angle = sill;

    # distance from points to center of sill
    Δ0 = p - Center

    # rotate point
    Δ =  InjectSills.rotate_point(Δ0, sill.RotMat.val)

    # sum of squares of distances (done as loop to avoid allocations)
    r = zero(_T)
    for i=1:N-1
        r += Δ[i]^2
    end
    r = sqrt(r)
    z = abs(Δ[N])

    # Function below cannot deal with zero
    if r==0; r=_T(1e-8); end
    if z==0; z=_T(1e-8); end

    # Compute displacement, using complex functions
    Ur, Uz = compute_penny_shaped_displacement_complex(r, z, ΔP, ν, E, R)

    if (Δ[N]<0); Uz = -Uz; end

    # Ur is the axisymmetric radial (outward) displacement; project onto the in-plane axes
    Displacement = N == 2 ? Vec2{_T}(Δ[1]/r*Ur, Uz) : Vec3{_T}(Δ[1]/r*Ur, Δ[2]/r*Ur, Uz)

    # rotate backwards
    Displacement_r = InjectSills.rotate_point(Displacement, sill.RotMat.val')


    return Displacement_r
end


function compute_penny_shaped_displacement(r, z, ΔP, ν, E, R)
    # Attempt to mimic the complex implementation of Sun (below) using real numbers
    error("This does not work correctly, use the complex implementation")
    # WRONG!!!
    R1 = sqrt(r^2 + (z - R)^2)
    R2 = sqrt(r^2 + (z + R)^2)

    # equation 7a:
    term1 = r * log((R2 + z + R) / (R1 + z - R))
    term2 = r / 2 * ((R - 3z - R2) / (R2 + z + R) + (R1 + 3z + R) / (R1 + z - R))
    term3 = (2 * z^2 * r) / (1 - 2ν) * (1 / (R2 * (R2 + z + R)) - 1 / (R1 * (R1 + z - R)))
    term4 = (2 * z * r) / (1 - 2ν) * (1 / R2 - 1 / R1)
    dU = ΔP * (1 + ν) * (1 - 2ν) / (2π * E) * (term1 - term2 - term3 + term4)

    # equation 7b:
    term5 = z * log((R2 + z + R) / (R1 + z - R))
    term6 = R2 - R1
    term7 = 1 / (2 * (1 - ν)) * (z * log((R2 + z + R) / (R1 + z - R)) - R * z * (1 / R2 + 1 / R1))
    dW = 2 * ΔP * (1 - ν^2) / (π * E) * (term5 - term6 - term7)

    Uz = dW  # vertical displacement should be corrected for z<0
    Ur = dU

    return Ur, Uz
end


# Complex division. Base divides Complex{Float32} in Float64, which GPUs without
# Float64 support (Metal) cannot compile.
cdiv(a, b) = a / b
cdiv(a, b::Complex{Float32}) = a * conj(b) / abs2(b)

"""
    compute_penny_shaped_displacement_complex(r, z, ΔP, ν, E, R)

Compute the displacement around a penny-shaped sill in a homogeneous elastic halfspace using the original implementation of Sun that uses complex numbers.
"""
function compute_penny_shaped_displacement_complex(r, z, ΔP, ν, E, R)
    imR = im*R
    R1  = sqrt(r^2 + (z - imR)^2);
    R2  = sqrt(r^2 + (z + imR)^2);
    L   = log(cdiv(R2+z+imR, R1+z-imR))

    # equation 7a:
    dU  = im*ΔP*(1+ν)*(1-2ν)/(2*E*π)*( r*L
            - r/2*(cdiv(imR-3z-R2, R2+z+imR)
            + cdiv(R1+3z+imR, R1+z-imR))
            - (2z^2 * r)/(1 -2ν)*(cdiv(1, R2*(R2+z+imR)) - cdiv(1, R1*(R1+z-imR)))
            + (2*z*r)/(1-2ν)*(cdiv(1, R2) - cdiv(1, R1)) );

    # equation 7b:
    dW  = 2*im*ΔP*(1-ν^2)/(pi*E)*( z*L
            - (R2-R1)
            - 1/(2*(1-ν))*( z*L - imR*z*(cdiv(1, R2) + cdiv(1, R1))) );

    Uz =  real(dW);  # vertical displacement should be corrected for z<0
    Ur =  real(dU);

    return Ur, Uz
end


"""
    inside(p::Point{2, _T}, sill::PennyShapedSill{2,_T})
checks if a 2D point `p` is inside the sill
"""
function inside(p::Point{2, _T}, sill::PennyShapedSill{2,_T}; rotate::Bool=true) where {_T}
    GeoParams.@unpack_val R,H, Center, RotMat = sill;

    # shift and optionally rotate point
    p_r = p - Center
    if rotate
        p_r = rotate_point(p_r, RotMat)
    end

    distance    = sqrt((p_r[1] / R)^2 + (p_r[2] / (H/2))^2)

    return distance <= 1
end



"""
    inside(p::Point{3, _T}, sill::PennyShapedSill{3,_T})
checks if a 3D point `p` is inside the sill
"""
function inside(p::Point{3, _T}, sill::PennyShapedSill{3,_T}; rotate::Bool=true) where {_T}
    GeoParams.@unpack_val R,H, Center, RotMat = sill;

    # shift and optionally rotate point
    p_r = p - Center
    if rotate
        p_r = rotate_point(p_r, RotMat)
    end

    distance = (p_r[1] / R)^2 + (p_r[2] / R)^2 + (p_r[3] / (H/2))^2

    return distance <= 1
end


"""
    update_abstractsill(s::PennyShapedSill; kwargs...) -> PennyShapedSill

Return a new sill identical to `s` but with the specified parameters updated.
Accepted keyword arguments: `Center`, `Angle`, `E`, `ν`, `R`, `H`, `ΔP`, `Q`.

When a single geometric parameter is changed, the remaining ones are kept
consistent using these defaults:
- only `R` → keep `H`, recompute `ΔP` and `Q`  (R+H branch)
- only `H` → keep `R`, recompute `ΔP` and `Q`  (R+H branch)
- only `ΔP` → keep `H`, recompute `R` and `Q`  (H+ΔP branch)
- only `Q` → keep `H`, recompute `R` and `ΔP`  (H+Q branch)

# Example
```julia
p2 = update_abstractsill(p, R = 3000.0m)
p3 = update_abstractsill(p, R = 3000.0m, H = 200.0m)
p4 = update_abstractsill(p, ΔP = 2e6Pa)
```
"""
function update_abstractsill(s::PennyShapedSill; kwargs...)
    reject_W(PennyShapedSill, "radius", get(kwargs, :W, nothing))
    check_keywords(PennyShapedSill, kwargs, (:Center, :Angle, :E, :ν, :R, :H, :ΔP, :Q))
    kw = Dict{Symbol,Any}(kwargs)

    has_R  = haskey(kw, :R)
    has_H  = haskey(kw, :H)
    has_ΔP = haskey(kw, :ΔP)
    has_Q  = haskey(kw, :Q)

    R  = get(kw, :R,  nothing)
    H  = get(kw, :H,  nothing)
    ΔP = get(kw, :ΔP, nothing)
    Q  = get(kw, :Q,  nothing)

    # Ensure a valid parameter pair reaches the constructor.
    if !has_R && !has_H && !has_ΔP && !has_Q
        R = UnitValue(s.R);  H = UnitValue(s.H)   # nothing changed: keep R+H
    elseif has_R && !has_H && !has_ΔP && !has_Q
        H = UnitValue(s.H)                          # only R: R+H branch
    elseif has_H && !has_R && !has_ΔP && !has_Q
        R = UnitValue(s.R)                          # only H: R+H branch
    elseif has_ΔP && !has_R && !has_H && !has_Q
        H = UnitValue(s.H)                          # only ΔP: H+ΔP branch
    elseif has_Q && !has_R && !has_H && !has_ΔP
        H = UnitValue(s.H)                          # only Q: H+Q branch
    # Otherwise user supplied a complete valid pair; pass as-is.
    end

    return PennyShapedSill(;
        Center = get(kw, :Center, UnitValue(s.Center)),
        Angle  = get(kw, :Angle,  UnitValue(s.Angle)),
        E      = get(kw, :E,      UnitValue(s.E)),
        ν      = get(kw, :ν,      UnitValue(s.ν)),
        R, H, ΔP, Q,
    )
end
