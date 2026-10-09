"""
    dike_polygon(sill::AbstractSill, nump::Integer=101)

Return a plotting polygon for supported kinematic dike-like sill types: the outline of the
sill in its own x–z plane through the center, rotated with the sill.
The polygon is returned as `[x, z]` in 2D and `[x, y, z]` in 3D, where `x`, `y` and `z` are vectors.
"""
function dike_polygon(sill::AbstractSill, nump::Integer=101)
    error("dike_polygon not implemented for type $(typeof(sill))")
end

# outline points (x[i], z[i]) in the x–z plane of the sill, relative to its center → model coordinates
function sill_to_model(sill::AbstractSill{2, _T}, x, z) where {_T}
    poly = [rotate_point(Point2{_T}(x[i], z[i]), sill.RotMat.val') + sill.Center.val for i in eachindex(x, z)]
    return [[pt[1] for pt in poly], [pt[2] for pt in poly]]
end

function sill_to_model(sill::AbstractSill{3, _T}, x, z) where {_T}
    poly = [rotate_point(Point3{_T}(x[i], zero(_T), z[i]), sill.RotMat.val') + sill.Center.val for i in eachindex(x, z)]
    return [[pt[1] for pt in poly], [pt[2] for pt in poly], [pt[3] for pt in poly]]
end

function dike_polygon(sill::Union{CylindricalDikeTopAccretion, CylindricalDikeTopAccretionFullModelAdvection}, nump::Integer=101)
    GeoParams.@unpack_val W, H = sill
    n = max(2, Int(nump))
    xx = collect(range(-W / 2, stop=W / 2, length=n))
    x = [xx; xx[end:-1:1]]
    z = [fill(H / 2, n); fill(-H / 2, n)]
    return sill_to_model(sill, x, z)
end

function dike_polygon(sill::EllipticalIntrusion{N, _T}, nump::Integer=101) where {N, _T}
    GeoParams.@unpack_val W, H = sill
    n = max(4, Int(nump))
    p = range(zero(_T), stop=2 * π, length=n)
    return sill_to_model(sill, cos.(p) .* (W / 2), -sin.(p) .* (H / 2))
end

function dike_polygon(sill::Union{PennyShapedSill, PlaneStrainSill}, nump::Integer=101)
    GeoParams.@unpack_val R, H = sill
    n = max(4, Int(nump))
    p = range(zero(R), stop=2 * π, length=n)
    return sill_to_model(sill, cos.(p) .* R, sin.(p) .* (H / 2))
end

function dike_polygon(sill::MogiSphere{2, _T}, nump::Integer=101) where {_T}
    GeoParams.@unpack_val r, Center = sill
    n = max(4, Int(nump))
    p = range(zero(_T), stop=2 * π, length=n)
    x = Center[1] .+ r .* cos.(p)
    z = Center[2] .+ r .* sin.(p)
    return [collect(x), collect(z)]
end

function dike_polygon(sill::McTigueSphere{2, _T}, nump::Integer=101) where {_T}
    GeoParams.@unpack_val r, Center = sill
    n = max(4, Int(nump))
    p = range(zero(_T), stop=2 * π, length=n)
    x = Center[1] .+ r .* cos.(p)
    z = Center[2] .+ r .* sin.(p)
    return [collect(x), collect(z)]
end
