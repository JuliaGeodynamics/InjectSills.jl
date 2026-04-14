
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

This generates a single new point that is within the sill
"""
function new_point_inside_sill(sill::AbstractSill{N,_T}) where {_T,N}
    GeoParams.@unpack_val W,H, Center, Angle = sill;

    isinside = false
    count = 1;
    coord = zero(Center)
    
    while !isinside && count<1000
        count += 1

        # generate a point that is within a square region that could be within the sill    
        coord = point_within_box(coord,W,H)

        # Rotate the point
        coord =  rotate_point(coord, sill.RotMat.val')   

        # Shift
        coord += Center

        # check if the point is within the sill
        isinside = inside(coord, sill)
        #isinside = true
    end
    

    return coord
end

"""
    pt = new_point_inside_sill(sill::AbstractSill{N,_T}, xvi, nx, ny; parts_per_cell = 10)

This generates a 3D Array with parts_per_cell-particles in each cell within the sill.
"""
function new_point_inside_sill(sill::AbstractSill{N,_T}, xvi, nx, ny; parts_per_cell = 10) where {_T,N}

    GeoParams.@unpack_val W,H, Center, Angle = sill

    # this layout kind-of mimics the layout of a CellArray
    parts2inject = fill(Point{N,_T}(NaN,NaN), nx, ny, parts_per_cell)
    dummy        = zero(Center)

    # iterate over cells
    for j in axes(parts2inject, 2), i in axes(parts2inject, 1)

        iscell_inside = rectangles_intercept(
            (xvi[1][i], xvi[2][j], xvi[1][i+1], xvi[2][j+1]),         # x1_min, y1_min, x1_max, y1_max
            (Center[1] - W/2, Center[2] - H/2, Center[1] + W/2, Center[2] + H/2)  # x2_min, y2_min, x2_max, y2_max
        )

        iscell_inside || continue

        # iterate over particles in the cell
        for k in axes(parts2inject, 3)
            # @show (i,j,k)

            isinside = false
            while !isinside
                # generate a point that is within a square region that could be within the sill    
                coord = point_within_box(dummy, W, H)
                # Rotate the point
                coord = rotate_point(coord, sill.RotMat.val')   
                # Shift
                coord += Center
                # check if the point is within the sill
                isinside = inside(coord, sill)
                if isinside
                    parts2inject[i, j, k] = coord
                end
            end

        end
    end

    return parts2inject
end

function rectangles_intercept(rect1, rect2)
    # Unpack the rectangles
    x1_min, y1_min, x1_max, y1_max = rect1
    x2_min, y2_min, x2_max, y2_max = rect2

    # Check if there is no overlap
    if x1_max < x2_min || x2_max < x1_min || y1_max < y2_min || y2_max < y1_min
        return false
    end

    # Otherwise, they intersect
    return true
end

point_within_box(p::Point{2, _T},W,H) where {_T} = Point2{_T}( 2*(rand()-0.5)*W, (rand()-0.5)*H )
point_within_box(p::Point{3, _T},W,H) where {_T} = Point3{_T}( 2*(rand()-0.5)*W, 2*(rand()-0.5)*W, (rand()-0.5)*H )


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
