using JustPIC
using JustPIC._2D
using InjectSills

# Threads is the default backend, 
# to run on a CUDA GPU load CUDA.jl (i.e. "using CUDA"), 
# and to run on an AMD GPU load AMDGPU.jl (i.e. "using AMDGPU")
const backend = JustPIC.CPUBackend # Options: CPUBackend, CUDABackend, AMDGPUBackend

using GLMakie

# 2D sill injection
function inject_sill!(particles,Dx,Dy, dx,dy, sill::AbstractSill{2,_T}) where {_T}

    # compute displacement field at the particles
    for I in eachindex(Dx.data)
        displacement = hostrock_displacement(sill, Point2{Float64}(px[I],py[I]))    
        Dx.data[I], Dy.data[I] = displacement[1], displacement[2]
    end

    fac = 1
    N   = 1
    if  sill.H.val>min(dx,dy)
        N = ceil(sill.H.val/min(dx,dy))
        fac = 1.0/N
    end
    @show fac


    for i=1:N
        # simple advection
        particles.coords[1].data .+= Dx.data*fac
        particles.coords[2].data .+= Dy.data*fac

        # we need to call `move_particles!` to update the particle data structure in case they moved to another cell
        move_particles!(particles, xvi, (Dx, Dy))
    end



    
    return nothing
end



#function main()
    # Initialize particles -------------------------------
    nxcell, max_xcell, min_xcell = 24, 30, 12
    n  = 256
    nx = ny = n-1
    Lx = Ly = 10000.0
    # nodal vertices
    xvi = xv, yv = range(-Lx/2, Lx/2, length=n), range(-Ly, 0, length=n)
    dxi = dx, dy = xv[2] - xv[1], yv[2] - yv[1]

    particles = init_particles(
        backend, nxcell, max_xcell, min_xcell, xvi...,
        buffer=1.0
    )

    # unpack particle coordinates
    px = particles.coords[1].data;
    py = particles.coords[2].data;

  
    # plot particles & their displacement
    #fig, ax, sc1 = scatter(px[:]/1e3, py[:]/1e3, markersize = 5, color=:black)
    
    # specify a sill 
    sill2D = PennyShapedSill(Center=Point2(0,-5000)*m, H=40.0m, W=200.0m, Angle=Vec1(0))

    # allocate particle arrays to hold the sill displacement field
    Dx, Dy = init_cell_arrays(particles, Val(2));
   
    # call the inject sill routine, which moves the particles
    inject_sill!(particles,Dx,Dy, dx,dy, sill2D) 

    inside_vec = zeros(Int32, size(px))
    for I in eachindex(px[:])
        inside_vec[I] = inside(Point2(px[I],py[I]), sill2D)
    end


    fig, ax, sc1 = scatter(px[:]/1e3, py[:]/1e3, color=Dy.data[:], markersize = 5)

    ind = findall(inside_vec.==1);
    scatter!(ax,px[ind]/1e3, py[ind]/1e3, markersize = 15)
    
    
   # inside(p::Point{2, _T}, sill::PennyShapedSill{2,_T})
    


    
    # plot ellipse
    # Number of points to create the ellipse
    N = 100
    theta = range(0, 2π, length=N)
     
    # Parametric equations for the ellipse
    x = (2*sill2D.W.val / 2) * cos.(theta) 
    y = (sill2D.H.val / 2) * sin.(theta)
    polygon =  Point2.(x,y)
    for I in eachindex(polygon)
        polygon[I] = InjectSills.rotate_point(polygon[I], sill2D.RotMat.val') + sill2D.Center.val
    end
    

    lines!(ax, polygon./1e3, color=:red)


    display(fig)

    # at this stage, we need to inject new particles within the sill, which is an ellipse region in 2D and an ellipsoide in 3D
    Ninject = 10000
    px_new,py_new = zeros(Ninject), zeros(Ninject)
    for i=1:Ninject
        px_new[i], py_new[i] = new_point_inside_sill(sill2D)
    end
    
  
    # plot newly added particles
    #scatter!(ax,px_new[:]/1e3, py_new[:]/1e3, markersize = 5, color=:green)
    
    display(fig)

#end


#main()