module InjectSillsJustPICExt

using InjectSills
using JustPIC, JustPIC._2D

"""
    inject_sill!(particles, Dx, Dy, xvi, sill::AbstractSill{2,_T})

Advect JustPIC `particles` by the displacement field induced by `sill`.
`Dx` and `Dy` are pre-allocated `CellArray`s (from `init_cell_arrays`).
`xvi` is the tuple of nodal-vertex ranges `(xv, yv)`.
"""
function InjectSills.inject_sill!(
    particles, Dx, Dy, xvi, sill::AbstractSill{2,_T}
) where {_T}
    dx = xvi[1][2] - xvi[1][1]
    dy = xvi[2][2] - xvi[2][1]

    px = particles.coords[1].data
    py = particles.coords[2].data

    # Compute displacement at every particle position
    for I in eachindex(Dx.data)
        d = hostrock_displacement(sill, Point2{_T}(px[I], py[I]))
        Dx.data[I] = d[1]
        Dy.data[I] = d[2]
    end

    # Sub-step if the sill is thicker than one cell
    N   = 1
    fac = 1.0
    if sill.H.val > min(dx, dy)
        N   = ceil(Int, sill.H.val / min(dx, dy)) * 2
        fac = 1.0 / N
    end

    for _ in 1:N
        particles.coords[1].data .+= Dx.data .* fac
        particles.coords[2].data .+= Dy.data .* fac
        move_particles!(particles, xvi, (Dx, Dy))
    end

    return nothing
end

end # module
