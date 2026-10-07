module InjectSillsJustPICExt

using InjectSills
using InjectSills: local_to_world
using JustPIC
using KernelAbstractions, Adapt

"""
    inject_sill!(particles, Dx, Dy, sill::AbstractSill{2}; fields=(), values=(), force_inject=false)
    inject_sill!(particles, Dx, Dy, Dz, sill::AbstractSill{3}; fields=(), values=(), force_inject=false)

Displace JustPIC `particles` by the host-rock displacement of `sill` and move them to
their new cells. `Dx`, `Dy` (, `Dz`) are `CellArray`s from `init_cell_arrays`; they
return the displacement of each particle. `fields` are further particle `CellArray`s
(e.g. phase, temperature) that move with the particles. The grid is taken from
`particles`. The displacements are computed with `hostrock_displacement!`.

With `force_inject=true`, the sill is then filled with new particles at the initial
density `particles.nxcell` per cell, and each of `fields` is set to the matching entry
of `values` on them.

All work runs on the backend of `particles` (CPU or GPU). On a GPU, the sill's float type
must be supported by the device, e.g. `Float32` on Metal.
"""
InjectSills.inject_sill!(particles, Dx, Dy, sill::AbstractSill{2}; kwargs...) =
    _inject_sill!(particles, (Dx, Dy), sill; kwargs...)

InjectSills.inject_sill!(particles, Dx, Dy, Dz, sill::AbstractSill{3}; kwargs...) =
    _inject_sill!(particles, (Dx, Dy, Dz), sill; kwargs...)

function _inject_sill!(particles, D::NTuple{N}, sill::AbstractSill{N, _T};
                       fields=(), values=(), force_inject=false) where {N, _T}
    force_inject && length(fields) != length(values) &&
        throw(ArgumentError("inject_sill!: fields and values must have the same length"))
    coords = map(c -> c.data, particles.coords)
    # empty particle slots have NaN coordinates
    hostrock_displacement!(map(d -> d.data, D), sill, coords; skipnan=true)

    for k in 1:N
        coords[k] .+= D[k].data
    end
    move_particles!(particles, (D..., fields...))

    if force_inject
        # new particles were not displaced: D = 0 on them
        force_injection!(particles, sill_particles(particles, sill), (fields..., D...), (values..., ntuple(_ -> zero(_T), N)...))
    end
    return nothing
end

# New particles inside `sill`, slot-aligned with `particles.index` as `force_injection!`
# expects, on the backend of `particles`. Each cell overlapping the sill draws `nxcell`
# uniform candidates and keeps those inside the sill, which reproduces the initial particle
# density, and places them in the cell's free slots.
function sill_particles(particles, sill::AbstractSill{N, _T}) where {N, _T}
    (; index, nxcell, xvi) = particles
    backend = get_backend(index.data)
    p_new = KernelAbstractions.allocate(backend, Point{N, _T}, size(index)..., cellnum(index))
    fill!(p_new, Point{N, _T}(ntuple(_ -> _T(NaN), Val(N))))

    # Cells overlapping the sill's bounding box; cell i spans xv[i]..xv[i+1]
    lo_s, hi_s = world_bounding_box(sill)
    xv = map(Array, xvi)
    cells = CartesianIndices(ntuple(Val(N)) do k
        max(searchsortedlast(xv[k], lo_s[k]), 1):min(searchsortedfirst(xv[k], hi_s[k]) - 1, length(xv[k]) - 1)
    end)
    isempty(cells) && return p_new

    # candidate positions within each cell, as fractions of the cell size
    U = adapt(backend, rand(_T, N, nxcell, size(cells)...))
    offset = first(cells) - oneunit(first(cells))
    dropped = KernelAbstractions.zeros(backend, Int, size(cells)...)
    sill_particles_kernel!(backend)(p_new, dropped, index, xvi, U, adapt(backend, sill), offset; ndrange = size(cells))
    synchronize(backend)
    n = sum(dropped)
    n == 0 || error("inject_sill!: $n new sill particles found no free slot in their cell; increase the particles' `max_xcell`")
    return p_new
end

@kernel function sill_particles_kernel!(p_new, dropped, index, xvi, U, sill, offset)
    J = @index(Global, Cartesian)
    dropped[J] = fill_cell!(p_new, index, xvi, U, sill, J, offset)
end

# Returns the number of candidates inside the sill that found no free slot.
function fill_cell!(p_new, index, xvi, U, sill::AbstractSill{N, _T}, J, offset) where {N, _T}
    I  = Tuple(J + offset)
    lo = Point{N, _T}(ntuple(k -> xvi[k][I[k]], Val(N)))
    hi = Point{N, _T}(ntuple(k -> xvi[k][I[k] + 1], Val(N)))
    ip = 0
    dropped = 0
    for c in axes(U, 2)
        p = Point{N, _T}(ntuple(k -> lo[k] + U[k, c, Tuple(J)...] * (hi[k] - lo[k]), Val(N)))
        inside(p, sill) || continue
        # next free slot (JustPIC.doskip is true for an empty slot)
        ip += 1
        while ip <= cellnum(index) && !JustPIC.doskip(index, ip, I...)
            ip += 1
        end
        if ip > cellnum(index)
            dropped += 1
        else
            p_new[I..., ip] = p
        end
    end
    return dropped
end

# World-frame axis-aligned box enclosing the (rotated) sill.
function world_bounding_box(sill::AbstractSill{N}) where {N}
    lo, hi = sill.BoundingBox[1].val, sill.BoundingBox[2].val
    corners = [local_to_world(sill, typeof(lo)(ntuple(k -> isodd(c >> (k - 1)) ? hi[k] : lo[k], Val(N))))
               for c in 0:(2^N - 1)]
    return reduce((a, b) -> min.(a, b), corners), reduce((a, b) -> max.(a, b), corners)
end

end # module
