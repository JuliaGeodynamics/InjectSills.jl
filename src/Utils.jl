
# few helper routines such as rotation matrixes 
using StaticArrays, GeometryBasics
export new_point_inside_sill

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



"""
    pt = new_point_inside_sill(sill::AbstractSill{N,_T})

Generates a single new point that is within the sill.
Samples uniformly from the unrotated bounding box in the sill frame, accepts the
first point that passes `inside(...; rotate=false)`, and maps it to world
coordinates with `local_to_world`. Throws if no point is accepted after 1000 tries.
"""
function new_point_inside_sill(sill::AbstractSill{N,_T}) where {_T,N}
    lower = sill.BoundingBox[1].val
    upper = sill.BoundingBox[2].val
    for _ in 1:1000
        q = random_point_in_bbox(lower, upper)
        inside(q, sill; rotate=false) && return local_to_world(sill, q)
    end
    error("new_point_inside_sill: no point inside $(nameof(typeof(sill))) after 1000 tries")
end

# Maps a point of the unrotated, centered sill to world coordinates; sources without
# `RotMat` (spheres) are rotation invariant.
local_to_world(sill::AbstractSill, q) =
    hasproperty(sill, :RotMat) ? sill.Center.val + rotate_point(q - sill.Center.val, sill.RotMat.val') : q

# Sample a uniformly random point inside an axis-aligned bounding box.
random_point_in_bbox(lower::Point{2,_T}, upper::Point{2,_T}) where {_T} =
    Point2{_T}(lower[1] + rand(_T) * (upper[1] - lower[1]),
               lower[2] + rand(_T) * (upper[2] - lower[2]))

random_point_in_bbox(lower::Point{3,_T}, upper::Point{3,_T}) where {_T} =
    Point3{_T}(lower[1] + rand(_T) * (upper[1] - lower[1]),
               lower[2] + rand(_T) * (upper[2] - lower[2]),
               lower[3] + rand(_T) * (upper[3] - lower[3]))


# Build an axis-aligned (unrotated) bounding box around a center point.
function unrotated_bounding_box(center::GeoUnit{Point{2, _T}, U}, hx::_T, hz::_T) where {_T, U}
    c = center.val
    lower = convert(GeoUnit, Point2{_T}(c[1] - hx, c[2] - hz) * center.unit)
    upper = convert(GeoUnit, Point2{_T}(c[1] + hx, c[2] + hz) * center.unit)
    return (lower, upper)
end

function unrotated_bounding_box(center::GeoUnit{Point{3, _T}, U}, hx::_T, hy::_T, hz::_T) where {_T, U}
    c = center.val
    lower = convert(GeoUnit, Point3{_T}(c[1] - hx, c[2] - hy, c[3] - hz) * center.unit)
    upper = convert(GeoUnit, Point3{_T}(c[1] + hx, c[2] + hy, c[3] + hz) * center.unit)
    return (lower, upper)
end


# Create a named tuple from a struct, which is useful for some of the dispatches in the sill constructor
to_nt(s) = NamedTuple{fieldnames(typeof(s))}(Tuple(getfield(s, f) for f in fieldnames(typeof(s))))
