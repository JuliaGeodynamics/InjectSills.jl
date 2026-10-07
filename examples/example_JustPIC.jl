# Requires InjectSills, JustPIC and GLMakie in the active environment.
using InjectSills
using JustPIC
using GLMakie

const backend = JustPIC.CPU
const HOST, MAGMA = 1.0, 2.0

# Active particles: x, y, then the values of each of `fields`.
function active(particles, fields...)
    ok = vec(particles.index.data)
    return particles.coords[1].data[ok], particles.coords[2].data[ok], map(f -> f.data[ok], fields)...
end

function main()
    # max_xcell needs headroom: cells next to the opened sill end up with ~1.5 nxcell particles
    nxcell, max_xcell, min_xcell = 24, 48, 12
    n  = 256
    Lx = Ly = 10000.0
    xv = range(-Lx/2, Lx/2, length=n)
    yv = range(-Ly, 0, length=n)
    xvi = (xv, yv)

    # Staggered velocity grids: vertex vector on the diagonal, extended cell centers off it.
    xc  = xv[1:end-1] .+ step(xv) / 2
    yc  = yv[1:end-1] .+ step(yv) / 2
    xce = range(xc[1] - step(xv), xc[end] + step(xv), length=length(xc) + 2)
    yce = range(yc[1] - step(yv), yc[end] + step(yv), length=length(yc) + 2)
    particles = init_particles(backend, nxcell, max_xcell, min_xcell, (xv, yce), (xce, yv))

    sill = PennyShapedSill(Center=Point2(0, -5000)*m, H=40.0m, W=2000.0m, Angle=Vec1(30))

    Dx, Dy, phase = init_cell_arrays(particles, Val(3))
    phase.data .= HOST
    x0, y0, _ = active(particles, phase)

    inject_sill!(particles, Dx, Dy, xvi, sill; fields=(phase,), values=(MAGMA,), force_inject=true)

    x1, y1, ph, ux, uy = active(particles, phase, Dx, Dy)
    host, magma = ph .== HOST, ph .== MAGMA
    println("particles before: $(length(x0)), after: $(length(x1)) (host $(count(host)), magma $(count(magma)))")
    println("host particles inside the sill: $(count(i -> inside(Point2(x1[i], y1[i]), sill), findall(host)))")
    println("magma particles outside the sill: $(count(i -> !inside(Point2(x1[i], y1[i]), sill), findall(magma)))")

    # Before / after close-up around the sill center, and the full sill after injection
    xp, zp = dike_polygon(sill, 200)
    fig = Figure(size=(1500, 550))
    win = (-300, 300, -5300, -4700)
    for (col, (title, x, y, c)) in enumerate((("before", x0, y0, fill(HOST, length(x0))),
                                              ("after", x1, y1, ph)))
        ax = Axis(fig[1, col]; title, aspect=DataAspect(), xlabel="x [m]", ylabel="z [m]", limits=win)
        sel = @. win[1] <= x <= win[2] && win[3] <= y <= win[4]
        scatter!(ax, x[sel], y[sel]; color=c[sel], colorrange=(HOST, MAGMA), colormap=[:gray70, :red], markersize=3)
        lines!(ax, xp, zp; color=:black)
    end
    ax = Axis(fig[1, 3]; title="host-rock displacement |u| [m]", aspect=DataAspect(), xlabel="x [km]", ylabel="z [km]")
    sel = findall(host)[1:20:end]
    sc  = scatter!(ax, x1[sel] ./ 1e3, y1[sel] ./ 1e3; color=hypot.(ux[sel], uy[sel]), colormap=:viridis, markersize=2)
    lines!(ax, xp ./ 1e3, zp ./ 1e3; color=:red)
    Colorbar(fig[1, 4], sc)

    save(joinpath(@__DIR__, "example_JustPIC.png"), fig)
    display(fig)
    return particles, phase
end

main()
