# example that shows how we can download a topography of a certain region 
# and compute the surface deformation cause by dike injection.
using GMT, GeophysicalModelGenerator, InjectSills



"""
    Topo_cart_disp = surface_deformation_sill(Topo_cart::CartData, sill::AbstractSill{3,_T}) where {_T}

Computes surface deformation caused by a sill and adds that as a field to `Topo_cart`. 
"""
function surface_deformation_sill(Topo_cart::CartData, sill::AbstractSill{3,_T}) where {_T}
    sz = size(Topo_cart)
    Dx = zeros(_T, sz)
    Dy = zeros(_T, sz)
    Dz = zeros(_T, sz)
        
    # compute displacement field at the particles
    for I in eachindex(Topo_cart.x.val)
        p = Point3{_T}(Topo_cart.x.val[I]*_T(1e3), Topo_cart.y.val[I]*_T(1e3), Topo_cart.z.val[I]*_T(1e3))
        displacement = hostrock_displacement(sill, p)    
        Dx[I], Dy[I], Dz[I] = displacement[1], displacement[2], displacement[3]
    end
    
    surface_deformation_m = (Dx,Dy,Dz)

    return addfield(Topo_cart, (; surface_deformation_m,))
end


# add a routine that creates a 3D VTK file of the sill.




# GMT command:
#grdimage("@earth_relief_01s", region=(-91.5,-90.9,-1,-0.6), coast=true, fmt=:png, show=true)

# download topography (of part of Galapagos):
Topo = import_topo([-91.2,-91.05,-0.9,-0.75], file="@earth_relief_03s");

# Create a cartesian grid from this
proj      = ProjectionPoint(; Lat=mean(Topo.lat.val[:]), Lon=mean(Topo.lon.val[:]))
Topo_cart = convert2CartData(Topo, proj)

# Create a sill. Note that the dimensions are in meters!
sill3D = PennyShapedSill(Center=Point3(0,0,-3000)*m, H=100.0m, W=2000.0m, Angle=Vec2(10,-20))

# Compute surface deformation caused by the sill
Topo_cart = surface_deformation_sill(Topo_cart, sill3D)

# 3D displacement field
XYZ             = xyz_grid(-5:.05:5,-5:.05:5,-10:.05:2)
Data            = CartData(XYZ)
displacement    = hostrock_displacement(sill3D, XYZ[1]*1000, XYZ[2]*1000, XYZ[3]*1000)
Data            = addfield(Data, (;displacement,))

# specify a sill
write_paraview(Topo_cart, "Topo_galapagos")
write_paraview(Data, "Displacement_galapagos")
