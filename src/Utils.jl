"""
    pt = new_point_inside_sill(sill::AbstractSill{N,_T}, xvi, nx, ny; parts_per_cell = 10)

This generates a 3D Array with parts_per_cell-particles in each cell within the sill.
"""
function new_point_inside_sill(sill::AbstractSill{N,_T}, xvi, nx, ny; parts_per_cell = 10) where {_T,N}

    GeoParams.@unpack_val W,H, Center, Angle = sill

    # this layout kind-of mimics the layout of a CellArray
    parts2inject = fill(Point{N,_T}(NaN,NaN), nx, ny, parts_per_cell)
    dummy        = zero(Center)

    for j in axes(parts2inject, 2)[1:end-1], i in axes(parts2inject, 1)[1:end-1]

        iscell_inside = any(
            inside(Point{N,_T}(xvi[1][ii], xvi[2][jj]), sill) for ii in i:i+1, jj in j:j+1
        )

        iscell_inside || continue

        for k in axes(parts2inject, 3)[1:end-1]
            # @show (i,j,k)

            isinside = false
            while !isinside
                # generate a point that is within a square region that could be within the sill    
                coord = point_within_box(dummy, W, H)
                # Rotate the point
                coord = rotate_point(coord, sill.RotMat.val')   
                # Shift
                coord += Center
                # check if the point is within the sill
                isinside = inside(coord, sill)
                if isinside
                    parts2inject[i, j, k] = coord
                end
            end

        end
    end

    return parts2inject
end
