"""
    dike_polygon(sill::AbstractSill, nump::Integer=101)

Return a 2D plotting polygon for supported kinematic dike-like sill types.
The polygon is returned as `[x, z]`, where `x` and `z` are vectors.
"""
function dike_polygon(sill::AbstractSill, nump::Integer=101)
    error("dike_polygon not implemented for type $(typeof(sill))")
end

function dike_polygon(sill::CylindricalDikeTopAccretion{N, _T}, nump::Integer=101) where {N, _T}
    GeoParams.@unpack_val W, H, Center = sill
    n = max(2, Int(nump))
    xx = collect(range(Center[1] - W / 2, stop=Center[1] + W / 2, length=n))
    x = [xx; xx[end:-1:1]]
    z = [fill(Center[N] + H, n); fill(Center[N] - H / 2, n)]
    return [x, z]
end

function dike_polygon(sill::CylindricalDikeTopAccretionFullModelAdvection{N, _T}, nump::Integer=101) where {N, _T}
    GeoParams.@unpack_val W, H, Center = sill
    n = max(2, Int(nump))
    xx = collect(range(Center[1] - W / 2, stop=Center[1] + W / 2, length=n))
    x = [xx; xx[end:-1:1]]
    z = [fill(Center[N] + H, n); fill(Center[N] - H / 2, n)]
    return [x, z]
end

function dike_polygon(sill::EllipticalIntrusion{N, _T}, nump::Integer=101) where {N, _T}
    GeoParams.@unpack_val W, H, Center = sill
    n = max(4, Int(nump))
    p = range(zero(_T), stop=2 * π, length=n)
    x = Center[1] .+ cos.(p) .* (W / 2)
    z = Center[N] .- sin.(p) .* (H / 2)
    return [collect(x), collect(z)]
end
