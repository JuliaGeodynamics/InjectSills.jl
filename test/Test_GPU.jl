using Test
using InjectSills, JustPIC

# GPU_BACKEND is set by runtests.jl from `--backend=...`
@eval using $(Symbol(GPU_BACKEND))
ArrayT    = getproperty(@__MODULE__, Dict("CUDA" => :CuArray, "AMDGPU" => :ROCArray, "Metal" => :MtlArray)[GPU_BACKEND])
JPBackend = getproperty(@__MODULE__, Dict("CUDA" => :CUDABackend, "AMDGPU" => :ROCBackend, "Metal" => :MetalBackend)[GPU_BACKEND])
floats    = GPU_BACKEND == "Metal" ? (Float32,) : (Float64, Float32)   # Metal has no Float64

function gpu_test_sills(::Type{F}) where {F}
    C2, C3 = Point2{F}(0, -5000)*m, Point3{F}(0, 0, -5000)*m
    A1, A2 = Vec1{F}(30)*NoUnits, Vec2{F}(30, 45)*NoUnits
    sills = AbstractSill[]
    for (C, A) in ((C2, A1), (C3, A2))
        push!(sills,
            PennyShapedSill(Center=C, Angle=A, R=F(1000)m, H=F(10)m, E=F(1.5e10)Pa, ν=F(0.3)*NoUnits),
            SquareDike(Center=C, Angle=A, W=F(2000)m, H=F(100)m),
            SquareDikeTopAccretion(Center=C, Angle=A, W=F(2000)m, H=F(100)m),
            CylindricalDikeTopAccretion(Center=C, Angle=A, W=F(2000)m, H=F(100)m),
            CylindricalDikeTopAccretionFullModelAdvection(Center=C, Angle=A, W=F(2000)m, H=F(100)m),
            EllipticalIntrusion(Center=C, Angle=A, W=F(2000)m, H=F(200)m),
            MogiSphere(Center=C, r=F(1500)m, ΔP=F(10e6)Pa, G=F(10e9)Pa, ν=F(0.25)*NoUnits),
            McTigueSphere(Center=C, r=F(1500)m, ΔP=F(10e6)Pa, G=F(10e9)Pa, ν=F(0.25)*NoUnits))
    end
    push!(sills, FiniteEllipsoidalCavity(Center=C3, ax=F(500)m, ay=F(300)m, az=F(2000)m,
        Angle=Vec{3,F}(10, 20, 45)*NoUnits, ΔP=F(1e6)Pa, mu=F(10e9)Pa, lambda=F(10e9)Pa))
    push!(sills, PlaneStrainSill(Center=C2, Angle=A1, R=F(1000)m, H=F(10)m, E=F(1.5e10)Pa, ν=F(0.3)*NoUnits))
    return sills
end

@testset "hostrock_displacement on $GPU_BACKEND, $F" for F in floats
    x = F[x for x in range(-3000, 3000, length=64), _ in 1:48]
    y = F[y for _ in 1:64, y in range(-2000, 2000, length=48)]
    for s in gpu_test_sills(F)
        N = length(s.Center.val)
        z = fill(F(s isa FiniteEllipsoidalCavity ? 0 : -4800), size(x))
        X = N == 2 ? (x, z) : (x, y, z)
        Dc = hostrock_displacement(s, X...)
        Dg = hostrock_displacement(s, map(ArrayT, X)...)
        @test all(D -> D isa ArrayT, Dg)
        @test all(k -> isapprox(Array(Dg[k]), Dc[k]; rtol = sqrt(eps(F))), 1:N)
    end
end

# Staggered velocity grids for a uniform vertex grid `xv`.
function gpu_staggered(xvs...)
    xcs = map(xv -> (xc = xv[1:end-1] .+ step(xv) / 2; range(xc[1] - step(xv), xc[end] + step(xv), length=length(xc) + 2)), xvs)
    return ntuple(d -> ntuple(e -> d == e ? xvs[e] : xcs[e], length(xvs)), length(xvs))
end

@testset "inject_sill! on $GPU_BACKEND, $(N)D" for N in (2, 3)
    F = first(floats)
    xvi = N == 2 ? (range(F(-5000), F(5000), length=65), range(F(-10000), F(0), length=65)) :
                   (range(F(-5000), F(5000), length=17), range(F(-5000), F(5000), length=17), range(F(-10000), F(0), length=17))
    particles = init_particles(JPBackend, 12, 36, 6, gpu_staggered(xvi...)...)
    sill = N == 2 ?
        PennyShapedSill(Center=Point2{F}(0, -5000)*m, H=F(40)m, R=F(2000)m, Angle=Vec1{F}(30)*NoUnits, E=F(1.5e10)Pa, ν=F(0.3)*NoUnits) :
        PennyShapedSill(Center=Point3{F}(0, 0, -5000)*m, H=F(400)m, R=F(2000)m, Angle=Vec2{F}(0, 0)*NoUnits, E=F(1.5e10)Pa, ν=F(0.3)*NoUnits)
    D..., phase = init_cell_arrays(particles, Val(N + 1))
    phase.data .= 1
    @test inject_sill!(particles, D..., sill; fields=(phase,), values=(F(2),), force_inject=true) === nothing

    ok = Array(vec(particles.index.data))
    x  = ntuple(k -> Array(particles.coords[k].data)[ok], N)
    d  = ntuple(k -> Array(D[k].data)[ok], N)
    ph = Array(phase.data)[ok]
    host, magma = findall(==(1), ph), findall(==(2), ph)
    @test !isempty(magma) && all(i -> inside(Point{N, F}(ntuple(k -> x[k][i], N)), sill), magma)
    # host particles sit at x0 + u(x0) with x0 = x - D
    @test all(host) do i
        u = hostrock_displacement(sill, Point{N, F}(ntuple(k -> x[k][i] - d[k][i], N)))
        isapprox(collect(u), collect(ntuple(k -> d[k][i], N)); rtol = sqrt(eps(F)), atol = 1e-3)
    end
end
