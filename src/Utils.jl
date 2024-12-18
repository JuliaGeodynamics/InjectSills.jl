
# few helper routines such as rotation matrixes 
using StaticArrays, GeometryBasics

function RotationMatrix(Angle::Vec{1, _T})  where {_T}
    sinDipAngle, cosDipAngle  = sincosd(Angle[1])
    return SMatrix{2,2}([cosDipAngle -sinDipAngle; sinDipAngle cosDipAngle])
end

function RotationMatrix(Angle::Vec{2, _T})  where {_T}
    sinDipAngle, cosDipAngle  = sincosd(-Angle[1])
    sinStrikeAngle, cosStrikeAngle  = sincosd(Angle[2])
    
    roty =  SMatrix{3,3}([cosDipAngle 0 sinDipAngle ; 0 1 0 ; -sinDipAngle 0  cosDipAngle]);       
    rotz =  SMatrix{3,3}([cosStrikeAngle -sinStrikeAngle 0 ; sinStrikeAngle cosStrikeAngle 0 ; 0 0 1])

    return roty*rotz
end

function rotate_point(p::Point{2, _T}, RotMat::SMatrix{2,2,_T,4}) where {_T}
    pt   = SVector{2,_T}(p.data...)
    pt_r = RotMat*pt
    return Point2{_T}(pt_r[1], pt_r[2])
end

function rotate_point(p::Vec{2,_T}, RotMat::SMatrix{2,2,_T,4}) where {_T}
    pt   = SVector{2,_T}(p.data...)
    pt_r = RotMat*pt
    return Vec2{_T}(pt_r[1], pt_r[2])
end

function rotate_point(p::Point{3, _T}, RotMat::SMatrix{3,3,_T,9}) where {_T}
    pt   = SVector{3,_T}(p.data...)
    pt_r = RotMat*pt
    return Point3{_T}(pt_r[1], pt_r[2], pt_r[3])
end

function rotate_point(p::Vec{3,_T}, RotMat::SMatrix{3,3,_T,9}) where {_T}
    pt   = SVector{3,_T}(p.data...)
    pt_r = RotMat*pt
    return Vec3{_T}(pt_r[1], pt_r[2], pt_r[3])
end

"""
    dX,dY,dZ = hostrock_displacement(sill::AbstractSill{3,_T}, X::AbstractArray{_T,N},Y::AbstractArray{_T,N},Z::AbstractArray{_T,N})

Creates a 3D displacement field caused by a magma-filles sill intrusion at the points `X,Y,Z`
"""
function  hostrock_displacement(sill::AbstractSill{3,_T}, X::AbstractArray{_T,N},Y::AbstractArray{_T,N},Z::AbstractArray{_T,N}) where {N,_T}
    Dx = zero(X)
    Dy = zero(X)
    Dz = zero(X)
    
    for I in CartesianIndices(X)
        p = Point3{_T}(X[I], Y[I], Z[I])
        Dx[I], Dy[I], Dz[I] = hostrock_displacement(sill, p)    
    end

    return Dx, Dy, Dz
end

"""
    dX,dZ = hostrock_displacement(sill::AbstractSill{2,_T}, X::AbstractArray{_T,N}, Z::AbstractArray{_T,N})

Creates a 2D displacement field caused by a magma-filles sill intrusion at the points `X,Z`
"""
function  hostrock_displacement(sill::AbstractSill{2,_T}, X::AbstractArray{_T,N},Z::AbstractArray{_T,N}) where {N,_T}
    Dx = zero(X)
    Dz = zero(X)
    
    for I in CartesianIndices(X)
        p = Point2{_T}(X[I], Z[I])
        Dx[I], Dz[I] = hostrock_displacement(sill, p)    
    end

    return Dx, Dz
end
