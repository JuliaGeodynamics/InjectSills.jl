module InjectSillsJustPICExt

using InjectSills
using InjectSills: local_to_world, random_point_in_bbox
using JustPIC

"""
    inject_sill!(particles, Dx, Dy, xvi, sill::AbstractSill{2}; fields=(), values=(), force_inject=false)
    inject_sill!(particles, Dx, Dy, Dz, xvi, sill::AbstractSill{3}; fields=(), values=(), force_inject=false)

Displace JustPIC `particles` by the host-rock displacement of `sill` and move them to
their new cells. `Dx`, `Dy` (, `Dz`) are `CellArray`s from `init_cell_arrays`; they
return the displacement of each particle. `fields` are further particle `CellArray`s
(e.g. phase, temperature) that move with the particles. `xvi` is not used; the grid
is taken from `particles`.

With `force_inject=true`, the sill is then filled with new particles at the initial
density `particles.nxcell` per cell, and each of `fields` is set to the matching entry
of `values` on them.
"""
InjectSills.inject_sill!(particles, Dx, Dy, xvi, sill::AbstractSill{2}; kwargs...) =
    _inject_sill!(particles, (Dx, Dy), xvi, sill; kwargs...)

InjectSills.inject_sill!(particles, Dx, Dy, Dz, xvi, sill::AbstractSill{3}; kwargs...) =
    _inject_sill!(particles, (Dx, Dy, Dz), xvi, sill; kwargs...)

function _inject_sill!(particles, D::NTuple{N}, xvi, sill::AbstractSill{N, _T};
                       fields=(), values=(), force_inject=false) where {N, _T}
    coords = map(c -> c.data, particles.coords)
    for I in eachindex(coords[1])
        p = Point{N, _T}(ntuple(k -> coords[k][I], Val(N)))
        d = isnan(p) ? zero(Vec{N, _T}) : hostrock_displacement(sill, p)
        for k in 1:N
            D[k].data[I] = d[k]
        end
    end

    # move_particles! drops particles that transiently overfill a cell, which grows quickly
    # with the per-call displacement; substeps keep it below 1/16 of a cell.
    nsub = max(1, ceil(Int, 16 * maximum(k -> maximum(abs, D[k].data), 1:N) / minimum(particles.di.vertex)))
    for _ in 1:nsub
        for k in 1:N
            coords[k] .+= D[k].data ./ nsub
        end
        move_particles!(particles, (D..., fields...))
    end

    if force_inject
        length(fields) == length(values) || throw(ArgumentError("fields and values must have the same length"))
        # new particles were not displaced: D = 0 on them
        force_injection!(particles, sill_particles(particles, sill), (fields..., D...), (values..., ntuple(_ -> zero(_T), N)...))
    end
    return nothing
end

# New particles inside `sill`, slot-aligned with `particles.index` as `force_injection!`
# expects. Each cell overlapping the sill draws `nxcell` uniform candidates and keeps those
# inside the sill, which reproduces the initial particle density, and places them in the
# cell's free slots.
function sill_particles(particles, sill::AbstractSill{N, _T}) where {N, _T}
    (; index, nxcell, xvi) = particles
    p_new = fill(Point{N, _T}(ntuple(_ -> _T(NaN), Val(N))), size(index)..., cellnum(index))
    lo_s, hi_s = world_bounding_box(sill)
    for I in CartesianIndices(index)
        lo = Point{N, _T}(ntuple(k -> xvi[k][I[k]], Val(N)))
        hi = Point{N, _T}(ntuple(k -> xvi[k][I[k] + 1], Val(N)))
        (all(lo .<= hi_s) && all(lo_s .<= hi)) || continue
        free = findall(!, index[I])
        n = 0
        for _ in 1:nxcell
            p = random_point_in_bbox(lo, hi)
            inside(p, sill) || continue
            n += 1
            n > length(free) && break
            p_new[I, free[n]] = p
        end
    end
    return p_new
end

# World-frame axis-aligned box enclosing the (rotated) sill.
function world_bounding_box(sill::AbstractSill{N}) where {N}
    lo, hi = sill.BoundingBox[1].val, sill.BoundingBox[2].val
    corners = [local_to_world(sill, typeof(lo)(ntuple(k -> isodd(c >> (k - 1)) ? hi[k] : lo[k], Val(N))))
               for c in 0:(2^N - 1)]
    return reduce((a, b) -> min.(a, b), corners), reduce((a, b) -> max.(a, b), corners)
end

end # module
