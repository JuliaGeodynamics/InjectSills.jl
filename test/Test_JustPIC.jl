using Test
using JustPIC
using GeoParams, InjectSills

# Staggered velocity grids for a uniform vertex grid `xv`.
function staggered(xvs...)
    xcs = map(xv -> (xc = xv[1:end-1] .+ step(xv) / 2; range(xc[1] - step(xv), xc[end] + step(xv), length=length(xc) + 2)), xvs)
    return ntuple(d -> ntuple(e -> d == e ? xvs[e] : xcs[e], length(xvs)), length(xvs))
end

# Every host particle must sit at x0 + u(x0), where x0 = x - D is its position before injection.
advected(sill, x, D) = hostrock_displacement(sill, Point(x .- D)) ≈ Vec(D)

@testset "2D" begin
    xvi = (range(-5000.0, 5000.0, length=33), range(-10000.0, 0.0, length=33))
    particles = init_particles(JustPIC.CPU, 12, 30, 6, staggered(xvi...)...)
    sill = PennyShapedSill(Center=Point2(0, -5000)*m, H=40.0m, W=2000.0m, Angle=Vec1(30))
    Dx, Dy, phase = init_cell_arrays(particles, Val(3))
    phase.data .= 1.0
    @test inject_sill!(particles, Dx, Dy, xvi, sill; fields=(phase,), values=(2.0,), force_inject=true) === nothing

    ok = vec(particles.index.data)
    px, py, ph = particles.coords[1].data[ok], particles.coords[2].data[ok], phase.data[ok]
    dx, dy = Dx.data[ok], Dy.data[ok]
    host, magma = findall(==(1.0), ph), findall(==(2.0), ph)
    @test any(i -> !iszero(dy[i]), host)
    @test all(i -> advected(sill, (px[i], py[i]), (dx[i], dy[i])), host)
    @test !isempty(magma) && all(i -> inside(Point2(px[i], py[i]), sill), magma)
    @test_throws "fields and values" inject_sill!(particles, Dx, Dy, xvi, sill; fields=(phase,), force_inject=true)
end

@testset "3D" begin
    xvi = (range(-5000.0, 5000.0, length=17), range(-5000.0, 5000.0, length=17),
           range(-10000.0, 0.0, length=17))
    particles = init_particles(JustPIC.CPU, 6, 8, 3, staggered(xvi...)...)
    sill = PennyShapedSill(Center=Point3(0, 0, -5000)*m, H=400.0m, W=2000.0m, Angle=Vec2(0, 0))
    Dx, Dy, Dz = init_cell_arrays(particles, Val(3))
    @test inject_sill!(particles, Dx, Dy, Dz, xvi, sill; force_inject=true) === nothing
    ok = vec(particles.index.data)
    @test any(ok)
    @test all(isfinite, Dz.data[ok])
    @test any(i -> inside(Point3(particles.coords[1].data[i], particles.coords[2].data[i], particles.coords[3].data[i]), sill), findall(ok))
end
