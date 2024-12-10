# Penny shaped sills in an elastic halfspace

# Penny shapes sills embedded in a homogeneous elastic halfspace 
using Adapt

export PennyShapedSill

"""
    PennyShapedSill{_T,N}

Holds information about a penny shaped sill in 2D or 3D.

"""
@kwdef struct PennyShapedSill{_T,N}
    Angle       ::  NTuple{N, _T} =   ntuple(i -> 0.0, N)
    Type        ::  Symbol =   :ElasticDike
    T           ::  _T =   950.0
    E           ::  _T =   1.5e10
    ν           ::  _T =   0.3
    ΔP          ::  _T =   1e6
    Q           ::  _T =   1000
    W           ::  _T =   (3*E*Q/(16*(1-ν^2)*ΔP))^(1.0/3.0)
    H           ::  _T =   8*(1-ν^2)*ΔP*W/(π*E)
    Center      ::  NTuple{N, _T} =   ntuple(i -> 20e3, N)
    Phase       ::  Int64 =   2
end
Adapt.@adapt_structure PennyShapedSill





