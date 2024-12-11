# Penny shapes sills embedded in a homogeneous elastic halfspace 
using Adapt, StaticArrays
import Base.show
import GeoParams: isdimensional


export PennyShapedSill, set_penny_shaped_sill, hostrock_displacement


"""
    PennyShapedSill{N,N1,_T}

Holds information about a penny shaped sill in 2D or 3D.

Parameters:
====
- `Center::Point{N, _T}`  - Center of sill
- `Angle::Vec{N1, _T}`       - Dip and strike angle of sill w.r.t. horizontal
- `E::_T`                                - Young's modulus
- `ν::_T`                                   - Poisson's ratio
- `ΔP::_T`                                  - Overpressure within sill
- `Q::_T`                                - Total injected volume of sill [m^3]
- `W::_T = (3*E*Q/(16*(1-ν^2)*ΔP))^(1.0/3.0)`     - Width of sill
- `H::_T = 8*(1-ν^2)*ΔP*W/(π*E)`                  - Maximum thickness of sill


Reference:
===
   Sun, R.J., 1969. Theoretical size of hydraulically induced horizontal fractures and 
        corresponding surface uplift in an idealized medium. J. Geophys. Res. 74, 5995–6011. 
        https://doi.org/10.1029/JB074i025p05995

"""
struct PennyShapedSill{N, _T, N1, U1, U2, U3, U4, U5} <: AbstractSill{N,_T}
    Center::GeoUnit{Point{N, _T},U1}  # m
    Angle::GeoUnit{Vec{N1, _T},U2}    # degrees
    E::GeoUnit{_T,U3}   # in Pa   
    ν::GeoUnit{_T,U4}   # []
    ΔP::GeoUnit{_T,U3}  # Pa
    Q::GeoUnit{_T,U5}   # m^3  
    W::GeoUnit{_T,U1}   # m  
    H::GeoUnit{_T,U1}   # m
end
PennyShapedSill(args...) = PennyShapedSill(convert.(GeoUnit, args)...)
Adapt.@adapt_structure PennyShapedSill

isdimensional(PennyShapedSill) = isdimensional(PennyShapedSill.E)

"""
    PennyShapedSill(; W=nothing,  Q=nothing, ΔP=nothing, H=nothing, E=1.5e10Pa, ν=0.3*NoUnits, Angle=Vec1(0.0)*Pas, Center=Point2(0.0)*m)

Defines parameters for a penny shaped sill in an elastic halfspace.    
You can give various combinations of parameters to define the sill:
- `ΔP` and `H` 
- `H` and `Q`
- `W` and `H`
- `W` and `Q` 
"""
function PennyShapedSill(; W=nothing,  Q=nothing, ΔP=nothing, H=nothing, E=1.5e10Pa, ν=0.3*NoUnits, Angle=Vec1(0.0)*Pas, Center=Point2(0.0)*m)
    @assert length(Center)==length(Angle)+1

    if isnothing(W) && isnothing(Q) && isnothing(ΔP) && isnothing(H)
        ΔP  =  1e6*Pa
        Q   =  1000.0*m^3    
        W   =  (3*E*Q/(16*(1-ν^2)*ΔP))^(1.0/3.0)
        H   =  8*(1-ν^2)*ΔP*W/(π*E)

    elseif isnothing(W) && isnothing(Q) && !isnothing(ΔP) && isnothing(H)
        Q   =  1000.0*m^3    
        W   =  (3*E*Q/(16*(1-ν^2)*ΔP))^(1.0/3.0)
        H   =  8*(1-ν^2)*ΔP*W/(π*E)

    elseif isnothing(W) && !isnothing(Q) && !isnothing(ΔP) && isnothing(H)
        W   =  (3*E*Q/(16*(1-ν^2)*ΔP))^(1.0/3.0)
        H   =  8*(1-ν^2)*ΔP*W/(π*E)

    elseif isnothing(W) && !isnothing(Q) && isnothing(ΔP) && isnothing(H)
        ΔP  =  1e6*Pa
        W   =  (3*E*Q/(16*(1-ν^2)*ΔP))^(1.0/3.0)
        H   =  8*(1-ν^2)*ΔP*W/(π*E)

    elseif !isnothing(W) && !isnothing(Q)
        ΔP = 3*E*Q/(16*(1-ν^2)*W^3)
        H  = (3 * Q) / (2 * π * W^2)

    elseif !isnothing(W) && !isnothing(H)
        ΔP = (π * E * H) / (8 * (1 - ν^2) * W)
        Q = (2 * π * H * W^2) / 3

    elseif !isnothing(H) && !isnothing(Q)
        W = sqrt(3 * Q / (2 * π * H))
        ΔP = (π * E * H) / (8 * (1 - ν^2) * W)

    elseif !isnothing(H) && !isnothing(ΔP)
        W = (π * E * H) / (8 * (1 - ν^2) * ΔP)
        Q = (π^3 * E^2 * H^3) / (96 * (1 - ν^2)^2 * ΔP^2)

    end

    if isnothing(W) || isnothing(Q) || isnothing(ΔP) || isnothing(H)
        error("you need to specify W and Q or ΔP and H or combinations")
    end

    return PennyShapedSill(Center, Angle, E, ν, ΔP, Q, W, H)
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
    println(io, "   Maximum sill width      : $(UnitValue(g.W)) ")  
 
    return nothing
end


"""
rotate a 2D point `p` by an angle `angle` (given in degrees)
"""
function rotate_point!(p::Union{Point{2, _T}, Vec{2,_T}}, angle::Vec{1, _T}) where {_T}
    
    α       = angle[1]
    RotMat  = SMatrix{2,2}([cosd(α) -sind(α); sind(α) cosd(α)]);    # 2D rotation matrix
    p       = RotMat * p
    return nothing
end

"""
    set_penny_shaped_sill(sill::PennyShapedSill{N,_T}, args...)

host rock displacement caused by opening of a penny shaped sill at point `p`
"""
function hostrock_displacement(sill::PennyShapedSill{N,_T}, p::Point{N, _T}) where {N,_T}

    @unpack_val ν,E,W,H, ΔP, Center, Angle = sill;

    # distance from points to center of sill
    Δ = p - Center

    # rotate point
    #rotate_point!(Δ, Angle)    # allocates

    # 
    r = zero(_T)
    for i=1:N
        r += Δ[i]^2
    end
    r = sqrt(r)
    z = Δ[N]

    if r==0; r=1e-3; end

    # Compute displacement, using complex functions
    imW = im*W
    R1  = sqrt(r^2. + (z - imW)^2);
    R2  = sqrt(r^2. + (z + imW)^2);

    # equation 7a:
    dU  = im*ΔP*(1+ν)*(1-2ν)/(2pi*E)*( r*log( (R2+z+imW)/(R1 +z- imW)) 
            - r/2*((imW-3z-R2)/(R2+z+imW) 
            + (R1+3z+imW)/(R1+z-imW)) 
            - (2z^2 * r)/(1 -2ν)*(1/(R2*(R2+z+imW)) -1/(R1*(R1+z-imW))) 
            + (2*z*r)/(1-2ν)*(1/R2 - 1/R1) );

    # equation 7b:
    dW  = 2*im*ΔP*(1-ν^2)/(pi*E)*( z*log( (R2+z+imW)/(R1+z-imW)) 
            - (R2-R1) 
            - 1/(2*(1-ν))*( z*log( (R2+z+imW)/(R1+z-imW)) - imW*z*(1/R2 + 1/R1)) );

    Uz =  real(dW);  # vertical displacement should be corrected for z<0
    Ur =  real(dU);
    if (p[N]<0); Uz = -Uz; end
    if (p[1]<0); Ur = -Ur; end

    Displacement  = Vec{N, _T}(Uz)
    if N==2
        Displacement = Vec2(Ur,Uz)
    elseif N==3
        x,y   = abs.(p[1:2]); 
        Displacement = Vec3(x/r*Ur,y/r*Ur,Uz)
    end

    # rotate backwards
    # to be implemented

    return Displacement
end

