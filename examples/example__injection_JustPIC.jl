using JustPIC, JustPIC._2D
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
        N = ceil(sill.H.val/min(dx,dy))*2
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
    #sill2D = PennyShapedSill(Center=Point2(0,-5000)*m, H=250.0m, W=1000.0m, Angle=Vec1(0))
    sill2D = PennyShapedSill(Center=Point2(0,-5000)*m, H=15.0m, W=1000.0m, Angle=Vec1(0))
    


    # allocate particle arrays to hold the sill displacement field
    Dx, Dy = init_cell_arrays(particles, Val(2));
   
    # call the inject sill routine, which moves the particles
    inject_sill!(particles,Dx,Dy, dx,dy, sill2D) 

    fig, ax, sc1 = scatter(px[:]/1e3, py[:]/1e3, color=Dy.data[:], markersize = 5)
    
    
    # plot ellipse
    # Number of points to create the ellipse
    N = 100
    theta = range(0, 2π, length=N)
     
    # Parametric equations for the ellipse
    x = (2*sill2D.W.val / 2) * cos.(theta) .+  sill2D.Center[1].val
    y = (sill2D.H.val / 2) * sin.(theta) .+  sill2D.Center[2].val
    lines!(ax, x/1e3, y/1e3, color=:red)

    display(fig)

    p_new = new_point_inside_sill(sill2D, xvi, nx, ny; parts_per_cell = 50);

    # scatter!(ax, p_new./1e3, markersize = 5, color=:green)
    # display(fig)

    #### PROTO ####
    phases, = init_cell_arrays(particles, Val(1))
    phases.data .= 1

    force_injection!(particles, p_new, (phases,), (2e0,))
    #force_injection!(particles, p_new)

    # plot newly added particles
    scatter(
        particles.coords[1].data[:]./1e3, 
        particles.coords[2].data[:]./1e3, 
        color = phases.data[:]./1e3,
        markersize = 5, 
    )
    
    # display(fig)

#end


#main()
