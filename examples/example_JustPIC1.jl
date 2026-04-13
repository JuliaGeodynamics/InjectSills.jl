function main(sill2D:: PennyShapedSill)
    # Initialize particles -------------------------------
    nxcell, max_xcell, min_xcell = 24, 30, 12
    n  = 256
    nx = ny = n-1
    Lx = Ly = 30000.0
    # nodal vertices
    xvi = xv, yv = range(-Lx/2, Lx/2, length=n), range(-Ly, 0, length=n)
    dxi = dx, dy = xv[2] - xv[1], yv[2] - yv[1]

    particles = init_particles(
        backend, nxcell, max_xcell, min_xcell, xvi...,
    )

    # unpack particle coordinates
    px = particles.coords[1].data;
    py = particles.coords[2].data;

    # plot particles
    # fig, ax, sc1 = scatter(px[:]/1e3, py[:]/1e3, markersize = 5, color=:black)

    # specify a sill 
    #sill2D = PennyShapedSill(Center=Point2(0,-15000)*m, H=100.0m, W=10000.0m, Angle=Vec1(20))

    # @code_warntype PennyShapedSill(Center=Point2(0,-15000)*m, H=100.0m, W=10000.0m, Angle=Vec1(20))
    # @edit PennyShapedSill(Center=Point2(0,-15000)*m, H=100.0m, W=10000.0m, Angle=Vec1(20))
    # @b PennyShapedSill(Center=$(Point2(0,-15000)*m), H=$(100.0m), W=$(10000.0m), Angle=$(Vec1(20)))

    # allocate particle arrays to hold the sill displacement field
    Dx, Dy = init_cell_arrays(particles, Val(2));

    # compute displacement field at the particles
    allocs = 0
    for I in eachindex(Dx.data)
        allocs += @allocated hostrock_displacement(sill2D, Point2{Float64}(px[I], py[I]))    
        displacement = hostrock_displacement(sill2D, Point2{Float64}(px[I], py[I]))    
        Dx.data[I], Dy.data[I] = displacement[1], displacement[2]
    end

    # simple advection
    particles.coords[1].data .+= Dx.data
    particles.coords[2].data .+= Dy.data

    # we need to call `move_particles!` to update the particle data structure in case they moved to another cell
    move_particles!(particles, xvi, ())
    allocs
end

# specify a sill 
sill2D = PennyShapedSill(Center=Point2(0,-15000)*m, H=100.0m, W=10000.0m, Angle=Vec1(20))

main(sill2D)