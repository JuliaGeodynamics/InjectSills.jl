using Test
using InjectSills, Makie

for sill in (PennyShapedSill(Angle=Vec1(30.0)*NoUnits), PlaneStrainSill(Angle=Vec1(30.0)*NoUnits))
    fig, ax = plot_sill(sill; scale=1e-3, color=:red)
    @test fig isa Figure
    # the outline is rotated about the sill center
    pts = only(ax.scene.plots)[1][][1:end-1]   # the last point closes the outline
    @test sum(pts) / length(pts) ≈ sill.Center.val * 1e-3 atol=1e-6
    @test maximum(p -> hypot((p - sill.Center.val * 1e-3)...), pts) ≈ sill.R.val * 1e-3 rtol=1e-3
end
