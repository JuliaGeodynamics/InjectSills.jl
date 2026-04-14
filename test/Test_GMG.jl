using Test
using GeophysicalModelGenerator
using GeoParams, InjectSills

# ---- build a CartData surface following the README pattern -----------------
# GMG surfaces are rank-3 arrays with nz=1.  Coordinates are in km;
# InjectSills uses metres — the extension converts internally.

nx, ny, nz = 101, 101, 1
x2D = [xi for xi in range(-50.0, 50.0, length=nx), _ in 1:ny, _ in 1:nz]  # km
y2D = [yi for _ in 1:nx, yi in range(-50.0, 50.0, length=ny), _ in 1:nz]
z2D = zeros(nx, ny, nz)   # flat surface at z = 0 km
surf = CartData(x2D, y2D, z2D, (Elevation = z2D,))

# ---- MogiSphere: centre at (0, 0, -5 km) -----------------------------------

src = MogiSphere(
    Center = Point3(0.0, 0.0, -5000.0)*m,
    r  = 1500.0m,
    ΔP = 10e6Pa,
    G  = 10e9Pa,
    ν  = 0.25*NoUnits,
)

# ---- array return form ------------------------------------------------------

Ux, Uy, Uz = surface_displacement(src, surf)

@test size(Ux) == (nx, ny, nz)
@test size(Uy) == (nx, ny, nz)
@test size(Uz) == (nx, ny, nz)

# At (0, 0) — directly above centre — vertical uplift matches 1D Mogi result;
# horizontal components vanish by axisymmetry.
cx, cy = div(nx, 2) + 1, div(ny, 2) + 1
@test Uz[cx, cy, 1] ≈ 0.10125  rtol=1e-4
@test abs(Ux[cx, cy, 1]) < 1e-6
@test abs(Uy[cx, cy, 1]) < 1e-6

# Symmetry: Ux anti-symmetric along x, Uy along y
@test Ux[cx+10, cy, 1] ≈ -Ux[cx-10, cy, 1]  rtol=1e-4
@test Uy[cx, cy+10, 1] ≈ -Uy[cx, cy-10, 1]  rtol=1e-4

# ---- add_fields return form -------------------------------------------------
# Displacement stored as (Ux, Uy, Uz) tuple under :Displacement_m

surf_out = surface_displacement(src, surf; add_fields = true)

@test surf_out isa CartData
@test haskey(surf_out.fields, :Displacement_m)
@test haskey(surf_out.fields, :Elevation)   # original field preserved

# Components: Displacement_m[1]=Ux, [2]=Uy, [3]=Uz
@test surf_out.fields.Displacement_m[3][cx, cy, 1] ≈ Uz[cx, cy, 1]        rtol=1e-6
@test surf_out.fields.Displacement_m[1][cx+10, cy, 1] ≈ Ux[cx+10, cy, 1]  rtol=1e-6

# ---- topographic surface (z ≠ 0) -------------------------------------------
# Confirms z is converted to metres, not ignored.

z_topo = fill(-0.5, nx, ny, nz)   # 500 m below sea level (in km)
surf_topo = CartData(x2D, y2D, z_topo, (Elevation = z_topo,))

Ux_t, Uy_t, Uz_t = surface_displacement(src, surf_topo)

# At (0, 0, −500 m) source-to-surface distance is 4500 m vs 5000 m at z=0;
# vertical uplift must be larger.
@test Uz_t[cx, cy, 1] > Uz[cx, cy, 1]
