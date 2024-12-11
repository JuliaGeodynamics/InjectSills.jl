# Penny shapes sills embedded in a homogeneous elastic halfspace 
using Adapt, GeoParams
import Base.show
import GeoParams: isdimensional, UnitValue


export PennyShapedSill, set_penny_shaped_sill


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
struct PennyShapedSill{N, N1, _T, U1, U2, U3, U4, U5} 
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
        H  = 8*(1-ν^2)*ΔP*W/(π*E)

    elseif !isnothing(W) && !isnothing(H)
        ΔP = (π * E * H) / (8 * (1 - ν^2) * W)
        Q = (2 * π * H * W^2) / 3

        
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
    println(io, "   Maximum sill radius     : $(UnitValue(g.W)) ")  
 
    return nothing
end
