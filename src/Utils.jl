
# few helper routines such as rotation matrixes 
using StaticArrays, GeometryBasics, KernelAbstractions, Adapt
export new_point_inside_sill, hostrock_displacement!

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
    dX, dY, dZ = hostrock_displacement(sill::AbstractSill{3,_T}, X::AbstractArray{_T}, Y::AbstractArray{_T}, Z::AbstractArray{_T})

Displacement field of `sill` at the points `X`, `Y`, `Z`, computed with [`hostrock_displacement!`](@ref).
"""
function hostrock_displacement(sill::AbstractSill{3,_T}, X::AbstractArray{_T,N}, Y::AbstractArray{_T,N}, Z::AbstractArray{_T,N}) where {N,_T}
    return hostrock_displacement!((similar(X), similar(X), similar(X)), sill, (X, Y, Z))
end

"""
    dX, dZ = hostrock_displacement(sill::AbstractSill{2,_T}, X::AbstractArray{_T}, Z::AbstractArray{_T})

Displacement field of `sill` at the points `X`, `Z`, computed with [`hostrock_displacement!`](@ref).
"""
function hostrock_displacement(sill::AbstractSill{2,_T}, X::AbstractArray{_T,N}, Z::AbstractArray{_T,N}) where {N,_T}
    return hostrock_displacement!((similar(X), similar(X)), sill, (X, Z))
end

"""
    D = hostrock_displacement!(D, sill::AbstractSill{N}, X; skipnan=false)

Write the displacement of `sill` at the points `(X[1][i], …, X[N][i])` to
`(D[1][i], …, D[N][i])` for every index `i`, and return `D`. `D` and `X` are `N`-tuples of
arrays with identical axes. The points are evaluated in parallel by a KernelAbstractions
kernel on the backend of `X[1]`. With `skipnan=true`, points with a `NaN` coordinate get
zero displacement.
"""
function hostrock_displacement!(D::NTuple{N, AbstractArray}, sill::AbstractSill{N}, X::NTuple{N, AbstractArray}; skipnan::Bool=false) where {N}
    ax = axes(X[1])
    all(A -> axes(A) == ax, (X..., D...)) ||
        throw(DimensionMismatch("hostrock_displacement!: all arrays must have axes $ax; got $(map(axes, (X..., D...)))"))
    check_points(sill, X)
    backend = get_backend(X[1])
    # moves array fields of the sill (FiniteEllipsoidalCavity sources) to the backend
    displacement_kernel!(backend)(D, adapt(backend, sill), X, skipnan; ndrange = length(X[1]))
    synchronize(backend)
    return D
end

# Point displacement inside the kernel, where nothing may throw (GPU kernels cannot build
# error messages); `check_points` validates the inputs on the host before the launch.
kernel_displacement(sill, p) = hostrock_displacement(sill, p)
check_points(sill, X) = nothing

@kernel function displacement_kernel!(D, sill, X, skipnan)
    i = @index(Global, Linear)
    displacement_at!(D, sill, X, i, skipnan)
end

function displacement_at!(D, sill::AbstractSill{N, _T}, X, i, skipnan) where {N, _T}
    p = Point{N, _T}(ntuple(k -> X[k][i], Val(N)))
    d = skipnan && isnan(p) ? zero(Vec{N, _T}) : kernel_displacement(sill, p)
    ntuple(k -> (D[k][i] = d[k]), Val(N))
    return nothing
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


# Argument validation shared by the sill constructors and `update_abstractsill`.
# Values may be plain numbers, Unitful quantities, or `GeoUnit`s.
plain_value(x) = x isa GeoUnit ? UnitValue(x) : x

function check_positive(T, name, x)
    v = plain_value(x)
    ustrip(v) > 0 || throw(ArgumentError("$T: `$name` must be positive; got $name = $v"))
    return nothing
end

# Isotropic linear elasticity requires -1 < ν ≤ 1/2. Solutions containing
# 1/(1-2ν) are singular at ν = 1/2, so those callers pass `incompressible=false`.
function check_poisson_ratio(T, ν; incompressible::Bool)
    v = ustrip(plain_value(ν))
    valid = incompressible ? -1 < v <= 0.5 : -1 < v < 0.5
    valid || throw(ArgumentError("$T: Poisson's ratio must satisfy -1 < ν $(incompressible ? "≤" : "<") 0.5; got ν = $v"))
    return nothing
end

function check_keywords(T, kwargs, valid)
    for k in keys(kwargs)
        k in valid || throw(ArgumentError("$T: unknown keyword `$k`; valid keywords are $(join(valid, ", "))"))
    end
    return nothing
end
