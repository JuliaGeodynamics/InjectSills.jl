using InjectSills, Test

@testset "Penny shaped sill" begin
    include("PennyShapedSill.jl")
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

@testset "Finite Ellipsodial Cavity" begin
    include("Test_FiniteEllisoidalCavity.jl")
end

@testset "Mogi and McTigue" begin
    include("Test_Mogi_McTigue.jl")
end

try
    using GeophysicalModelGenerator
    @testset "GeophysicalModelGenerator extension" begin
        include("Test_GMG.jl")
    end
catch e
    @info "Skipping GMG extension tests" exception=e
end
