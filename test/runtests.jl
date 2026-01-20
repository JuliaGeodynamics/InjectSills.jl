using InjectSills, Test

@testset "Penny shaped sill" begin
    include("PennyShapedSill.jl")
end

@testset "Finite Ellipsodial Cavity" begin
    include("Test_FiniteEllisoidalCavity.jl")
end

@testset "Mogi and McTigue" begin
    include("Test_Mogi_McTigue.jl")
end