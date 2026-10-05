using Test
using JustPIC
using GeoParams, InjectSills

# Staggered velocity grids for a uniform vertex grid `xv`.
function staggered(xvs...)
    xcs = map(xv -> (xc = xv[1:end-1] .+ step(xv) / 2; range(xc[1] - step(xv), xc[end] + step(xv), length=length(xc) + 2)), xvs)
    return ntuple(d -> ntuple(e -> d == e ? xvs[e] : xcs[e], length(xvs)), length(xvs))
end

@testset "2D" begin
    xvi = (range(-5000.0, 5000.0, length=33), range(-10000.0, 0.0, length=33))
    particles = init_particles(JustPIC.CPU, 12, 18, 6, staggered(xvi...)...)
    sill = PennyShapedSill(Center=Point2(0, -5000)*m, H=40.0m, W=2000.0m, Angle=Vec1(30))
    Dx, Dy = init_cell_arrays(particles, Val(2))
    @test inject_sill!(particles, Dx, Dy, xvi, sill) === nothing

    px, py = particles.coords[1].data, particles.coords[2].data
    ok = .!isnan.(px) .& .!isnan.(py)
    @test any(ok)
    @test all(isfinite, Dx.data[ok]) && all(isfinite, Dy.data[ok])
    @test any(!iszero, Dy.data[ok])
    I = findfirst(ok)
    d = hostrock_displacement(sill, Point2(px[I], py[I]))
    @test (Dx.data[I], Dy.data[I]) == (d[1], d[2])
end

@testset "3D" begin
    xvi = (range(-5000.0, 5000.0, length=17), range(-5000.0, 5000.0, length=17),
           range(-10000.0, 0.0, length=17))
    particles = init_particles(JustPIC.CPU, 6, 8, 3, staggered(xvi...)...)
    sill = PennyShapedSill(Center=Point3(0, 0, -5000)*m, H=400.0m, W=2000.0m, Angle=Vec2(0, 0))
    Dx, Dy, Dz = init_cell_arrays(particles, Val(3))
    @test inject_sill!(particles, Dx, Dy, Dz, xvi, sill) === nothing
    ok = .!isnan.(particles.coords[1].data)
    @test any(ok)
    @test all(isfinite, Dz.data[ok])
end
