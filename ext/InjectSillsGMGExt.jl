module InjectSillsGMGExt

using InjectSills
using GeophysicalModelGenerator

"""
    Ux, Uy, Uz = surface_displacement(sill::AbstractSill{3}, cart::CartData)
    cart_out   = surface_displacement(sill::AbstractSill{3}, cart::CartData; add_fields=true)

Compute the surface displacement induced by `sill` at every point of the `CartData` surface
`cart`. The sill must be dimensional with lengths in m; otherwise an `ArgumentError` is thrown.
`CartData` coordinates are in km and are converted to meters internally before calling
`hostrock_displacement`.

Returns `(Ux, Uy, Uz)` as arrays in meters by default.

If `add_fields = true`, the displacement arrays are added to the `CartData` as named fields
`(:Ux, :Uy, :Uz)` in meters and the updated `CartData` is returned.

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
    # All length fields of a sill share the unit of `Center`.
    if !(InjectSills.isdimensional(sill.Center) && sill.Center.unit == m)
        got = InjectSills.isdimensional(sill.Center) ? "lengths in $(sill.Center.unit)" : "a nondimensionalized sill"
        throw(ArgumentError("surface_displacement: the sill must have lengths in m (CartData coordinates in km are converted to m); got $got"))
    end

    X = map(c -> _T.(1e3 .* c), (cart.x.val, cart.y.val, cart.z.val))   # km → m
    Ux, Uy, Uz = hostrock_displacement!(ntuple(_ -> similar(X[1]), 3), sill, X)

    if add_fields
        Displacement_m = (Ux, Uy, Uz)
        return addfield(cart, (; Displacement_m))
    else
        return Ux, Uy, Uz
    end
end

end # module InjectSillsGMGExt
