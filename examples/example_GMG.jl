# Example: download topography of a region and compute surface deformation
# caused by a sill/dike injection using InjectSills + GeophysicalModelGenerator.
using GMT, GeophysicalModelGenerator, InjectSills

# GMT command to preview the region (optional):
# grdimage("@earth_relief_01s", region=(-91.5,-90.9,-1,-0.6), coast=true, fmt=:png, show=true)

# Download topography (part of the Galapagos):
Topo = import_topo([-91.2, -91.05, -0.9, -0.75], file="@earth_relief_03s")

# Convert to a Cartesian grid
proj      = ProjectionPoint(; Lat=mean(Topo.lat.val[:]), Lon=mean(Topo.lon.val[:]))
Topo_cart = convert2CartData(Topo, proj)

# Define a sill (dimensions in metres)
sill3D = PennyShapedSill(Center=Point3(0, 0, -3000)*m, H=100.0m, W=2000.0m, Angle=Vec2(10, -20))

# Compute surface deformation and add as a field to the CartData
# Displacement_m = (Ux, Uy, Uz) tuple in metres
Topo_cart = surface_displacement(sill3D, Topo_cart; add_fields=true)

# 3D displacement field on a volume grid
XYZ          = xyz_grid(-5:0.05:5, -5:0.05:5, -10:0.05:2)
Data         = CartData(XYZ)
displacement = hostrock_displacement(sill3D, XYZ[1]*1000, XYZ[2]*1000, XYZ[3]*1000)
Data         = addfield(Data, (; displacement))

# Write to ParaView
write_paraview(Topo_cart, "Topo_galapagos")
write_paraview(Data, "Displacement_galapagos")
