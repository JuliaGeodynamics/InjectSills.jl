module InjectSillsGMGExt

using InjectSills
using GeophysicalModelGenerator

"""
    Ux, Uy, Uz = surface_displacement(sill::AbstractSill{3}, cart::CartData)
    cart_out   = surface_displacement(sill::AbstractSill{3}, cart::CartData; add_fields=true)

Compute the surface displacement induced by `sill` at every point of the `CartData` surface
`cart`. The sill must be in dimensional (SI) units. `CartData` coordinates are in km and are
converted to metres internally before calling `hostrock_displacement`.

Returns `(Ux, Uy, Uz)` as arrays in metres by default.

If `add_fields = true`, the displacement arrays are added to the `CartData` as named fields
`(:Ux, :Uy, :Uz)` in metres and the updated `CartData` is returned.

# Example
```julia
using InjectSills, GeophysicalModelGenerator

# flat surface at sea level
nx, ny = 201, 201
x2D = [xi for xi in range(-100.0, 100.0, length=nx), _ in 1:ny]  # km
y2D = [yi for _ in 1:nx, yi in range(-100.0, 100.0, length=ny)]
z2D = zeros(nx, ny)
surf = CartData(x2D, y2D, z2D, (Elevation=z2D,))

src = MogiSphere(
    Center = Point3(0.0, 0.0, -5000.0)*m,
    r = 1500.0m, ΔP = 10e6Pa, G = 10e9Pa, ν = 0.25*NoUnits,
)

# Option 1 — arrays
Ux, Uy, Uz = surface_displacement(src, surf)

# Option 2 — add as CartData fields
surf_out = surface_displacement(src, surf; add_fields=true)
```
"""
function InjectSills.surface_displacement(
    sill::AbstractSill{3, _T},
    cart::CartData;
    add_fields::Bool = false,
) where {_T}
    x_km = cart.x.val   # arrays in km
    y_km = cart.y.val
    z_km = cart.z.val

    Ux = similar(x_km, Float64)
    Uy = similar(x_km, Float64)
    Uz = similar(x_km, Float64)

    for I in eachindex(x_km)
        # convert km → m for the displacement calculation
        p = Point3{Float64}(
            Float64(x_km[I]) * 1e3,
            Float64(y_km[I]) * 1e3,
            Float64(z_km[I]) * 1e3,
        )
        d = hostrock_displacement(sill, p)
        Ux[I] = d[1]
        Uy[I] = d[2]
        Uz[I] = d[3]
    end

    if add_fields
        Displacement_m = (Ux, Uy, Uz)
        return addfield(cart, (; Displacement_m))
    else
        return Ux, Uy, Uz
    end
end

end # module InjectSillsGMGExt
