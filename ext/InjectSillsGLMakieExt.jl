module InjectSillsGLMakieExt

using InjectSills
using GLMakie

"""
    fig, ax = plot_sill(sill::PennyShapedSill{2}; scale=1.0, kwargs...)

Plot the outline of a 2-D penny-shaped sill.
`scale` rescales coordinates (e.g. `scale=1e-3` for km).
Extra `kwargs` are forwarded to `lines!`.
"""
function InjectSills.plot_sill(
    sill::PennyShapedSill{2,_T}; scale=1.0, kwargs...
) where {_T}
    N     = 200
    theta = range(0, 2π; length=N)
    x     = sill.W.val .* cos.(theta)
    y     = (sill.H.val / 2) .* sin.(theta)

    outline = Point2{_T}.(x, y)
    for I in eachindex(outline)
        outline[I] = InjectSills.rotate_point(outline[I], sill.RotMat.val') + sill.Center.val
    end

    fig = Figure()
    ax  = Axis(fig[1, 1]; aspect=DataAspect())
    lines!(ax, outline .* scale; kwargs...)
    return fig, ax
end

end # module
