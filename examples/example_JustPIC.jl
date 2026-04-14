import Pkg

function ensure_example_env!()
    # Keep JustPIC optional for InjectSills itself by using a local env for this example.
    env_dir = joinpath(@__DIR__, ".justpic_example_env")
    Pkg.activate(env_dir)

    pkgs = Set(keys(Pkg.project().dependencies))
    if !("InjectSills" in pkgs)
        Pkg.develop(path=joinpath(@__DIR__, ".."))
    end
    if !("JustPIC" in pkgs)
        Pkg.add(name="JustPIC")
    end

    Pkg.instantiate()
end

ensure_example_env!()

using InjectSills
using JustPIC
using JustPIC._2D

const backend = JustPIC.CPUBackend

function maybe_plot(px, py, Dy, sill2D)
    try
        @eval begin
            using GLMakie
        end

        fig, ax, _ = scatter(px[:]/1e3, py[:]/1e3, color=Dy.data[:], markersize=5)

        inside_vec = zeros(Int32, size(px))
        for I in eachindex(px)
            inside_vec[I] = inside(Point2(px[I], py[I]), sill2D)
        end
        ind = findall(inside_vec .== 1)
        scatter!(ax, px[ind]/1e3, py[ind]/1e3, markersize=15)

        n = 100
        theta = range(0, 2π, length=n)
        x = sill2D.W.val .* cos.(theta)
        y = (sill2D.H.val / 2) .* sin.(theta)
        polygon = Point2.(x, y)
        for I in eachindex(polygon)
            polygon[I] = InjectSills.rotate_point(polygon[I], sill2D.RotMat.val') + sill2D.Center.val
        end
        lines!(ax, polygon ./ 1e3, color=:red)

        display(fig)
    catch err
        @warn "Skipping plotting (GLMakie unavailable or no display)." err
    end
end

function main()
    nxcell, max_xcell, min_xcell = 24, 30, 12
    n = 256
    Lx = Ly = 10000.0
    xv = range(-Lx/2, Lx/2, length=n)
    yv = range(-Ly, 0, length=n)
    xvi = (xv, yv)

    particles = init_particles(backend, nxcell, max_xcell, min_xcell, xvi...)

    sill2D = PennyShapedSill(Center=Point2(0, -5000)*m, H=40.0m, W=200.0m, Angle=Vec1(0))

    Dx, Dy = init_cell_arrays(particles, Val(2))
    inject_sill!(particles, Dx, Dy, xvi, sill2D)

    px = particles.coords[1].data
    py = particles.coords[2].data

    println("JustPIC example ran successfully.")
    println("Particles: $(length(px))")
    println("Mean Dx: $(sum(Dx.data) / length(Dx.data))")
    println("Mean Dy: $(sum(Dy.data) / length(Dy.data))")

    maybe_plot(px, py, Dy, sill2D)
    return nothing
end

main()