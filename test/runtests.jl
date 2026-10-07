# `--backend=CUDA|AMDGPU|Metal` adds GPU tests (Buildkite); the default runs the CPU tests only.
backend_arg = filter(startswith("--backend="), ARGS)
GPU_BACKEND = isempty(backend_arg) ? "CPU" : split(only(backend_arg), "=")[2]
GPU_BACKEND in ("CPU", "CUDA", "AMDGPU", "Metal") ||
    error("unknown backend $(repr(GPU_BACKEND)); use --backend=CPU|CUDA|AMDGPU|Metal")

# GPU packages are not test dependencies, so the CPU tests never install them. The GPU package
# is added before any package is loaded, so that every loaded version comes from one resolve.
if GPU_BACKEND != "CPU"
    using Pkg
    Pkg.add(GPU_BACKEND)
end

using InjectSills, Test
using GeophysicalModelGenerator

@testset "Penny shaped sill" begin
    include("PennyShapedSill.jl")
end

@testset "Plane-strain sill" begin
    include("Test_PlaneStrainSill.jl")
end

@testset "Square dike sill" begin
    include("Test_SquareDike.jl")
end

@testset "Square dike top-accretion sill" begin
    include("Test_SquareDikeTopAccretion.jl")
end

@testset "Cylindrical dike top-accretion sill" begin
    include("Test_CylindricalDikeTopAccretion.jl")
end

@testset "Elliptical intrusion sill" begin
    include("Test_EllipticalIntrusion.jl")
end

@testset "Dike polygons" begin
    include("Test_DikePolygons.jl")
end

@testset "new_point_inside_sill" begin
    include("Test_NewPointInsideSill.jl")
end

@testset "Finite Ellipsodial Cavity" begin
    include("Test_FiniteEllisoidalCavity.jl")
end

@testset "Mogi and McTigue" begin
    include("Test_Mogi_McTigue.jl")
end

@testset "Array displacement" begin
    include("Test_ArrayDisplacement.jl")
end

@testset "Input validation" begin
    include("Test_InputValidation.jl")
end

@testset "GeophysicalModelGenerator extension" begin
    include("Test_GMG.jl")
end

@testset "JustPIC extension" begin
    include("Test_JustPIC.jl")
end

@testset "Makie extension" begin
    include("Test_Makie.jl")
end

if GPU_BACKEND != "CPU"
    @testset "GPU ($GPU_BACKEND)" begin
        include("Test_GPU.jl")
    end
end
