# Reference journal articles:
# 1)
# Nikkhoo, M., Rivalta, E. (2023):
# Surface deformations and gravity changes caused by pressurized finite 
# ellipsoidal cavities. Geophysical Journal International, doi: 
# 10.1093/gji/ggac351
#
# 2)
# Nikkhoo, M., Rivalta, E. (2022):
# Analytical solutions for gravity changes caused by triaxial volumetric
# sources. Geophysical Research Letters, doi: 10.1029/2021GL095442
#
# 3)
# Nikkhoo, M., Walter, T. R., Lundgren, P. R., Prats-Iraola, P. (2017):
# Compound dislocation models (CDMs) for volcano deformation analyses.
# Geophysical Journal International, doi: 10.1093/gji/ggw427
#
# 4)
# Eshelby, J. D. (1957):
# The determination of the elastic field of an ellipsoidal inclusion, and
# related problems.
# Proceedings of the royal society of London. Series A. Mathematical and
# physical sciences. 241 (1226), 376-396. doi: 10.1098/rspa.1957.0133
#
# 5)
# Carlson, B. C. (1995):
# Numerical computation of real or complex elliptic integrals.
# Numer. Algor., 10(1), 1326. doi: 10.1007/BF02198293
#
# 6)
# Klein, P. P. (2012):
# On the ellipsoid and plane intersection equation.
# Applied Mathematics, 3 (11), 1634-1640. doi: 10.4236/am.2012.311226

# Copyright (c) 2022 Mehdi Nikkhoo
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the
# following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
# NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.

# I appreciate any comments or bug reports.

# Mehdi Nikkhoo
# Created: 2021.7.24
# Last modified: 2022.10.24
#
# Section 2.1, Physics of Earthquakes and Volcanoes
# Department 2, Geophysics
# Helmholtz Centre Potsdam
# German Research Centre for Geosciences (GFZ)
#
# email:
# mehdi.nikkhoo@gfz-potsdam.de
# mehdi.nikkhoo@gmail.com
#
# website:
# http://www.volcanodeformation.com
#
# Translated to Julia by Arne Spang (2026.01.19)
# arne.spang@uni-bayreuth.de

using LinearAlgebra

"""
    fECM(X, Y, X0, Y0, depth, omegaX, omegaY, omegaZ, ax, ay, az, p, mu, lambda, DepthRef; Nmax=5e3, Cr=14)

Calculate the surface displacements caused by a uniformly-pressurized 
finite ellipsoidal cavity in a uniform elastic half-space.

fECM: finite Ellipsoidal Cavity Model
pCDM: point Compound Dislocation Model
PTD: Point Tensile Dislocation
EFCS: Earth-Fixed Coordinate System

# Arguments
- `X`, `Y`: Horizontal coordinates of calculation points in EFCS (East, North, Up).
  Can be scalars, vectors, or matrices. X and Y must have the same size.
- `X0`, `Y0`: Horizontal coordinates (in EFCS) of the center of the cavity.
- `depth`: The depth to the "reference point" (center or top of cavity).
  Positive number with same unit as X0 and Y0.
- `omegaX`, `omegaY`, `omegaZ`: Clockwise rotation angles about X, Y, Z axes (degrees).
- `ax`, `ay`, `az`: Semi-axes of the finite ECM along X, Y, Z axes before rotation.
- `p`: Pressure on the cavity walls (same unit as Lamé constants).
- `mu`, `lambda`: Lamé constants.
- `DepthRef`: Either 'C' or 'T' - reference point for depth.
  'C' = depth to center, 'T' = depth to top.

# Optional Arguments
- `Nmax`: Maximum total number of allowed point CDMs (default: 5e3)
- `Cr`: Grid spacing parameter (default: 14)

# Returns
- `ue`, `un`, `uv`: Displacement vector components in EFCS
- `dV`: Volume change of the finite ECM
- `DV`: Potency of the finite ECM
- `Ns`: Total number of point CDMs within the cavity

# Example
```julia
X, Y = 0:50:7000, zeros(length(0:50:7000))
X0, Y0 = 0.0, 0.0
depth = 2550.0
omegaX, omegaY, omegaZ = 11.0, -7.0, 120.0
ax, ay, az = 450.0, 600.0, 525.0
p = 1e6
mu, lambda = 1e9, 1e9
ue, un, uv, dV, DV, Ns = fECM(X, Y, X0, Y0, depth, omegaX, omegaY, omegaZ, 
                              ax, ay, az, p, mu, lambda, "C"; Nmax=4000, Cr=12)
```
"""
function fECM(X, Y, X0, Y0, depth, omegaX, omegaY, omegaZ,
              ax, ay, az, p, mu, lambda, DepthRef; 
              Nmax=5e3, Cr=14)
    
    # Get top and bottom of ellipsoid
    Zt, Zb = EllTopBot(depth, omegaX, omegaY, omegaZ, ax, ay, az)
    
    # Adjust depth based on reference point
    if uppercase(DepthRef) == "T"
        aV = (Zt - Zb) / 2
        depth = depth + aV  # Note that dC = depth + aV!
        Ztop = Zt - aV
    elseif uppercase(DepthRef) == "C"
        Ztop = Zt
    else
        error("Undefined input argument for 'DepthRef'! Must be 'C' or 'T'")
    end
    
    if Ztop >= 0
        error("Input error: the cavity is too shallow!")
    end
    
    # Flatten X and Y to vectors (handles scalars, vectors, and matrices)
    X_vec = X isa AbstractArray ? vec(X) : [X]
    Y_vec = Y isa AbstractArray ? vec(Y) : [Y]
    
    # Calculate material properties
    nu = lambda / (lambda + mu) / 2  # Poisson's ratio
    K = lambda + 2 * mu / 3  # Bulk modulus
    
    # Threshold for stable shape tensor calculation
    r0 = 1e-12
    ax = max(ax, r0)
    ay = max(ay, r0)
    az = max(az, r0)
    
    # Sort semi-axes in descending order
    ai = [ax, ay, az]
    perm = sortperm(ai, rev=true)
    ai_sorted = ai[perm]
    
    # Calculate shape tensor
    S = ShapeTensorECM(ai_sorted[1], ai_sorted[2], ai_sorted[3], nu)
    Sm = [S[1]-1  S[2]    S[3];
          S[4]    S[5]-1  S[6];
          S[7]    S[8]    S[9]-1]  # Shape tensor
    
    # Transformation strain
    eT = -inv(Sm) * p * ones(3) / (3 * K)
    
    # For uniformly-pressurized ellipsoids the eT elements have the same sign!
    eT[sign.(eT) .!= sign(p)] .= 0
    
    V = 4/3 * π * ax * ay * az  # Ellipsoid volume
    
    # Calculate actual volume change
    dV = (sum(eT) - p/K) * V
    
    # Calculate potency
    DV = dV + p * V / K  # Also DV = sum(M)/3/K could be used!
    
    # Reorder eT according to original axes order
    tmp = zeros(3)
    tmp[perm] = V * eT
    DVx = tmp[1]
    DVy = tmp[2]
    DVz = tmp[3]
    
    # Generate adaptive grid
    Xs, Ys, Zs, Ws = AdaptiveGrid(X0, Y0, depth, omegaX, omegaY, omegaZ, 
                                   ax, ay, az, Cr, Nmax)
    Ns = length(Xs)
    
    # Scale by weights
    DVx_scaled = DVx * Ws
    DVy_scaled = DVy * Ws
    DVz_scaled = DVz * Ws
    
    # Calculate displacements
    ue, un, uv = pCDM_Vec(X_vec, Y_vec, Xs, Ys, -Zs, omegaX, omegaY, omegaZ,
                          DVx_scaled, DVy_scaled, DVz_scaled, nu)
    
    # Reshape outputs to match input shape
    original_size = size(X)
    ue = reshape(ue, original_size)
    un = reshape(un, original_size)
    uv = reshape(uv, original_size)
    
    return ue, un, uv, dV, DV, Ns
end


"""
    pCDM_Vec(X, Y, X0, Y0, depth, omegaX, omegaY, omegaZ, DVx, DVy, DVz, nu)

Vectorized version of the pCDM (point Compound Dislocation Model) function.
Calculates surface displacements from compound dislocation sources.

# Arguments
- `X`, `Y`: Observation point coordinates (vectors)
- `X0`, `Y0`, `depth`: Source location coordinates (can be vectors for multiple sources)
- `omegaX`, `omegaY`, `omegaZ`: Rotation angles in degrees
- `DVx`, `DVy`, `DVz`: Potency components (can be vectors for multiple sources)
- `nu`: Poisson's ratio

# Returns
- `ue`, `un`, `uv`: Displacement components (East, North, Vertical)
"""
function pCDM_Vec(X, Y, X0, Y0, depth, omegaX, omegaY, omegaZ,
                  DVx, DVy, DVz, nu)
    
    # Ensure X and Y are vectors
    X = X isa AbstractArray ? vec(X) : [X]
    Y = Y isa AbstractArray ? vec(Y) : [Y]
    
    # Build rotation matrices
    Rx = [1  0            0;
          0  cosd(omegaX)  sind(omegaX);
          0 -sind(omegaX)  cosd(omegaX)]
    
    Ry = [cosd(omegaY)  0  -sind(omegaY);
          0             1   0;
          sind(omegaY)  0   cosd(omegaY)]
    
    Rz = [cosd(omegaZ)   sind(omegaZ)  0;
         -sind(omegaZ)   cosd(omegaZ)  0;
          0              0             1]
    
    R = Rz * Ry * Rx
    
    # Calculate strike and dip for first PTD (X-axis direction)
    Vstrike1 = [-R[2,1], R[1,1], 0]
    Vstrike1 = Vstrike1 / norm(Vstrike1)
    strike1 = atan(Vstrike1[1], Vstrike1[2]) * 180 / π
    if isnan(strike1)
        strike1 = 0.0
    end
    dip1 = acosd(R[3,1])
    
    # Calculate strike and dip for second PTD (Y-axis direction)
    Vstrike2 = [-R[2,2], R[1,2], 0]
    Vstrike2 = Vstrike2 / norm(Vstrike2)
    strike2 = atan(Vstrike2[1], Vstrike2[2]) * 180 / π
    if isnan(strike2)
        strike2 = 0.0
    end
    dip2 = acosd(R[3,2])
    
    # Calculate strike and dip for third PTD (Z-axis direction)
    Vstrike3 = [-R[2,3], R[1,3], 0]
    Vstrike3 = Vstrike3 / norm(Vstrike3)
    strike3 = atan(Vstrike3[1], Vstrike3[2]) * 180 / π
    if isnan(strike3)
        strike3 = 0.0
    end
    dip3 = acosd(R[3,3])
    
    # Calculate contribution of the first PTD
    ue1, un1, uv1 = PTD_disp_Surf_Vec(X, Y, X0, Y0, depth, strike1, dip1, DVx, nu)
    
    # Calculate contribution of the second PTD
    ue2, un2, uv2 = PTD_disp_Surf_Vec(X, Y, X0, Y0, depth, strike2, dip2, DVy, nu)
    
    # Calculate contribution of the third PTD
    ue3, un3, uv3 = PTD_disp_Surf_Vec(X, Y, X0, Y0, depth, strike3, dip3, DVz, nu)
    
    # Sum contributions
    ue = ue1 + ue2 + ue3
    un = un1 + un2 + un3
    uv = uv1 + uv2 + uv3
    
    return ue, un, uv
end

"""
Helper function to safely convert scalars or arrays to vectors
"""
vec_safe(x) = x isa AbstractArray ? vec(x) : [x]


"""
    PTD_disp_Surf_Vec(X, Y, X0, Y0, depth, strike, dip, DV, nu)

Vectorized version of PTD_disp_Surf function in the pCDM function.
Calculates surface displacement from Point Tensile Dislocation sources.

# Arguments
- `X`, `Y`: Observation point coordinates (vectors)
- `X0`, `Y0`, `depth`: Source location(s) - can be scalars or vectors
- `strike`, `dip`: Fault orientation parameters (degrees)
- `DV`: Potency - can be scalar or vector (same size as X0, Y0, depth)
- `nu`: Poisson's ratio

# Returns
- `ue`, `un`, `uv`: Displacement components (East, North, Vertical)
"""
function PTD_disp_Surf_Vec(X, Y, X0, Y0, depth, strike, dip, DV, nu)
    
    # Ensure X and Y are column vectors
    X = X isa AbstractArray ? vec(X) : [X]
    Y = Y isa AbstractArray ? vec(Y) : [Y]
    
    # X0, Y0, depth and DV must be row vectors of the same size
    # Use reshape to ensure they are proper 1×N matrices (not just transposed 1D arrays)
    X0 = X0 isa AbstractArray ? reshape(vec(X0), 1, :) : reshape([X0], 1, 1)
    Y0 = Y0 isa AbstractArray ? reshape(vec(Y0), 1, :) : reshape([Y0], 1, 1)
    depth = depth isa AbstractArray ? reshape(vec(depth), 1, :) : reshape([depth], 1, 1)
    DV = DV isa AbstractArray ? reshape(vec(DV), 1, :) : reshape([DV], 1, 1)
    
    # Create matrices using broadcasting (more efficient than repmat)
    xM = X .- X0  # Outer subtraction creates matrix
    yM = Y .- Y0
    d = ones(length(X), 1) * depth  # Column vector * row vector = matrix
    DV_mat = ones(length(X), 1) * DV
    
    beta = strike - 90
    x = xM .* cosd(beta) .- yM .* sind(beta)
    y = xM .* sind(beta) .+ yM .* cosd(beta)
    
    r = sqrt.(x.^2 .+ y.^2 .+ d.^2)
    q = y .* sind(dip) .- d .* cosd(dip)
    
    # Calculate integral terms
    I1 = (1 - 2*nu) .* y .* (1 ./ r ./ (r .+ d).^2 .- 
         x.^2 .* (3 .* r .+ d) ./ r.^3 ./ (r .+ d).^3)
    I2 = (1 - 2*nu) .* x .* (1 ./ r ./ (r .+ d).^2 .- 
         y.^2 .* (3 .* r .+ d) ./ r.^3 ./ (r .+ d).^3)
    I3 = (1 - 2*nu) .* x ./ r.^3 .- I2
    I5 = (1 - 2*nu) .* (1 ./ r ./ (r .+ d) .- 
         x.^2 .* (2 .* r .+ d) ./ r.^3 ./ (r .+ d).^2)
    
    # Note: For a PTD M0 = DV*mu!
    ue = DV_mat ./ (2 * π) .* (3 .* x .* q.^2 ./ r.^5 .- I3 .* sind(dip)^2)
    un = DV_mat ./ (2 * π) .* (3 .* y .* q.^2 ./ r.^5 .- I1 .* sind(dip)^2)
    uv = DV_mat ./ (2 * π) .* (3 .* d .* q.^2 ./ r.^5 .- I5 .* sind(dip)^2)
    
    # Sum over all sources
    ue0 = sum(ue, dims=2)
    un0 = sum(un, dims=2)
    uv = sum(uv, dims=2)
    
    # Rotate back to EFCS
    ue = ue0 .* cosd(beta) .+ un0 .* sind(beta)
    un = -ue0 .* sind(beta) .+ un0 .* cosd(beta)
    
    return vec(ue), vec(un), vec(uv)
end


"""
    ShapeTensorECM(a1, a2, a3, nu)

Calculate the Eshelby (1957) shape tensor components for an ellipsoid.

# Arguments
- `a1`, `a2`, `a3`: Semi-axes of ellipsoid (a1 ≥ a2 ≥ a3)
- `nu`: Poisson's ratio

# Returns
- `S`: 9-element vector containing shape tensor components [S1111, S1122, S1133, S2211, S2222, S2233, S3311, S3322, S3333]
"""
function ShapeTensorECM(a1, a2, a3, nu)
    
    if a1 == 0 && a2 == 0 && a3 == 0
        return zeros(9)
    end
    
    # Calculate Ik and Iij terms for triaxial, oblate and prolate ellipsoids
    if a1 > a2 && a2 > a3 && a3 > 0
        # General case: triaxial ellipsoid
        sin_theta = sqrt(1 - a3^2 / a1^2)
        k = sqrt((a1^2 - a2^2) / (a1^2 - a3^2))
        
        # Calculate Legendre's incomplete elliptic integrals using Carlson (1995) method
        tol = 1e-16
        c = 1 / sin_theta^2
        F = RF(c - 1, c - k^2, c, tol)
        E = F - k^2 / 3 * RD(c - 1, c - k^2, c, tol)
        
        I1 = 4 * π * a1 * a2 * a3 / (a1^2 - a2^2) / sqrt(a1^2 - a3^2) * (F - E)
        I3 = 4 * π * a1 * a2 * a3 / (a2^2 - a3^2) / sqrt(a1^2 - a3^2) *
             (a2 * sqrt(a1^2 - a3^2) / a1 / a3 - E)
        I2 = 4 * π - I1 - I3
        
        I12 = (I2 - I1) / (a1^2 - a2^2)
        I13 = (I3 - I1) / (a1^2 - a3^2)
        I11 = (4 * π / a1^2 - I12 - I13) / 3
        
        I23 = (I3 - I2) / (a2^2 - a3^2)
        I21 = I12
        I22 = (4 * π / a2^2 - I23 - I21) / 3
        
        I31 = I13
        I32 = I23
        I33 = (4 * π / a3^2 - I31 - I32) / 3
        
    elseif a1 == a2 && a2 > a3 && a3 > 0
        # Special case-1: Oblate ellipsoid
        I1 = 2 * π * a1 * a2 * a3 / (a1^2 - a3^2)^1.5 * 
             (acos(a3 / a1) - a3 / a1 * sqrt(1 - a3^2 / a1^2))
        I2 = I1
        I3 = 4 * π - 2 * I1
        
        I13 = (I3 - I1) / (a1^2 - a3^2)
        I11 = π / a1^2 - I13 / 4
        I12 = I11
        
        I23 = I13
        I22 = π / a2^2 - I23 / 4
        I21 = I12
        
        I31 = I13
        I32 = I23
        I33 = (4 * π / a3^2 - 2 * I31) / 3
        
    elseif a1 > a2 && a2 == a3 && a3 > 0
        # Special case-2: Prolate ellipsoid
        I2 = 2 * π * a1 * a2 * a3 / (a1^2 - a3^2)^1.5 * 
             (a1 / a3 * sqrt(a1^2 / a3^2 - 1) - acosh(a1 / a3))
        I3 = I2
        I1 = 4 * π - 2 * I2
        
        I12 = (I2 - I1) / (a1^2 - a2^2)
        I13 = I12
        I11 = (4 * π / a1^2 - 2 * I12) / 3
        
        I21 = I12
        I22 = π / a2^2 - I21 / 4
        I23 = I22
        
        I32 = I23
        I31 = I13
        I33 = (4 * π / a3^2 - I31 - I32) / 3
    end
    
    # Calculate the shape-tensor components
    if a1 == a2 && a2 == a3
        # Special case-3: Sphere
        S1111 = (7 - 5*nu) / 15 / (1 - nu)
        S1122 = (5*nu - 1) / 15 / (1 - nu)
        S1133 = (5*nu - 1) / 15 / (1 - nu)
        S2211 = (5*nu - 1) / 15 / (1 - nu)
        S2222 = (7 - 5*nu) / 15 / (1 - nu)
        S2233 = (5*nu - 1) / 15 / (1 - nu)
        S3311 = (5*nu - 1) / 15 / (1 - nu)
        S3322 = (5*nu - 1) / 15 / (1 - nu)
        S3333 = (7 - 5*nu) / 15 / (1 - nu)
    else
        # General triaxial, oblate and prolate ellipsoids
        S1111 = 3 / 8 / π / (1 - nu) * a1^2 * I11 + (1 - 2*nu) / 8 / π / (1 - nu) * I1
        S1122 = 1 / 8 / π / (1 - nu) * a2^2 * I12 - (1 - 2*nu) / 8 / π / (1 - nu) * I1
        S1133 = 1 / 8 / π / (1 - nu) * a3^2 * I13 - (1 - 2*nu) / 8 / π / (1 - nu) * I1
        S2211 = 1 / 8 / π / (1 - nu) * a1^2 * I21 - (1 - 2*nu) / 8 / π / (1 - nu) * I2
        S2222 = 3 / 8 / π / (1 - nu) * a2^2 * I22 + (1 - 2*nu) / 8 / π / (1 - nu) * I2
        S2233 = 1 / 8 / π / (1 - nu) * a3^2 * I23 - (1 - 2*nu) / 8 / π / (1 - nu) * I2
        S3311 = 1 / 8 / π / (1 - nu) * a1^2 * I31 - (1 - 2*nu) / 8 / π / (1 - nu) * I3
        S3322 = 1 / 8 / π / (1 - nu) * a2^2 * I32 - (1 - 2*nu) / 8 / π / (1 - nu) * I3
        S3333 = 3 / 8 / π / (1 - nu) * a3^2 * I33 + (1 - 2*nu) / 8 / π / (1 - nu) * I3
    end
    
    S = [S1111, S1122, S1133, S2211, S2222, S2233, S3311, S3322, S3333]
    
    return S
end


"""
    RF(x, y, z, r)

Calculate the RF term in the Carlson (1995) method for calculating elliptic integrals.

Carlson, B.C. (1995). Numerical computation of real or complex elliptic integrals.
Numerical Algorithms, 10, 13. doi:10.1007/BF02198293
"""
function RF(x, y, z, r)
    
    if any([x, y, z] .< 0)
        error("x, y and z values must be positive!")
    elseif count([x, y, z] .!= 0) < 2
        error("At most one of the x, y and z values can be zero!")
    end
    
    xm = x
    ym = y
    zm = z
    A0 = (x + y + z) / 3
    Q = maximum([abs(A0 - x), abs(A0 - y), abs(A0 - z)]) / (3 * r)^(1/6)
    n = 0
    Am = A0
    
    while abs(Am) <= Q / (4^n)
        lambdam = sqrt(xm * ym) + sqrt(xm * zm) + sqrt(ym * zm)
        Am = (Am + lambdam) / 4
        xm = (xm + lambdam) / 4
        ym = (ym + lambdam) / 4
        zm = (zm + lambdam) / 4
        n = n + 1
    end
    
    X = (A0 - x) / 4^n / Am
    Y = (A0 - y) / 4^n / Am
    Z = -X - Y
    E2 = X * Y - Z^2
    E3 = X * Y * Z
    rf = (1 - E2/10 + E3/14 + E2^2/24 - 3*E2*E3/44) / sqrt(Am)
    
    return rf
end


"""
    RD(x, y, z, r)

Calculate the RD term in the Carlson (1995) method for calculating elliptic integrals.
"""
function RD(x, y, z, r)
    
    if z == 0
        error("z value must be nonzero!")
    elseif x == 0 && y == 0
        error("At most one of the x and y values can be zero!")
    end
    
    xm = x
    ym = y
    zm = z
    A0 = (x + y + 3*z) / 5
    Q = maximum([abs(A0 - x), abs(A0 - y), abs(A0 - z)]) / (r / 4)^(1/6)
    n = 0
    Am = A0
    S = 0.0
    
    while abs(Am) <= Q / (4^n)
        lambdam = sqrt(xm * ym) + sqrt(xm * zm) + sqrt(ym * zm)
        S = S + (1 / 4^n) / sqrt(zm) / (zm + lambdam)
        Am = (Am + lambdam) / 4
        xm = (xm + lambdam) / 4
        ym = (ym + lambdam) / 4
        zm = (zm + lambdam) / 4
        n = n + 1
    end
    
    X = (A0 - x) / 4^n / Am
    Y = (A0 - y) / 4^n / Am
    Z = -(X + Y) / 3
    E2 = X * Y - 6 * Z^2
    E3 = (3 * X * Y - 8 * Z^2) * Z
    E4 = 3 * (X * Y - Z^2) * Z^2
    E5 = X * Y * Z^3
    rd = (1 - 3*E2/14 + E3/6 + 9*E2^2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26) / 
         4^n / Am^1.5 + 3 * S
    
    return rd
end


"""
    EllTopBot(depth, omegaX, omegaY, omegaZ, ax, ay, az)

Calculate the Z Cartesian coordinates of the shallowest (top) and deepest (bottom) 
points of an ellipsoid, plus their spherical coordinates.

The ellipsoid equation before applying rotations is: x²/ax² + y²/ay² + z²/az² = 1

# Returns
- `Ztop`, `Zbot`: Z coordinates of top and bottom points
- `Th`, `La`: Spherical coordinates (Theta and Lambda in degrees)
"""
function EllTopBot(depth, omegaX, omegaY, omegaZ, ax, ay, az)
    
    # Build rotation matrices
    Rx = [1  0            0;
          0  cosd(omegaX)  sind(omegaX);
          0 -sind(omegaX)  cosd(omegaX)]
    
    Ry = [cosd(omegaY)  0  -sind(omegaY);
          0             1   0;
          sind(omegaY)  0   cosd(omegaY)]
    
    Rz = [cosd(omegaZ)   sind(omegaZ)  0;
         -sind(omegaZ)   cosd(omegaZ)  0;
          0              0             1]
    
    R = Rz * Ry * Rx
    
    # Calculate spherical coordinates
    La = atan(R[3,2] * ay, R[3,1] * ax)
    Th = atan(R[3,1]^2 * ax^2 + R[3,2]^2 * ay^2,
              R[3,3] * az * sqrt(R[3,1]^2 * ax^2 + R[3,2]^2 * ay^2))
    
    if La < 0
        La = La + 2 * π
    end
    
    if Th < 0
        Th = Th + π
    elseif Th > π
        Th = Th - π
    end
    
    # Convert to degrees
    La = La * 180 / π
    Th = Th * 180 / π
    
    # Calculate Z coordinates
    Z = R[3,1] * ax * sind(Th) * cosd(La) +
        R[3,2] * ay * sind(Th) * sind(La) + 
        R[3,3] * az * cosd(Th)
    
    Z_vals = [Z, -Z] .- depth
    Ztop = maximum(Z_vals)
    Zbot = minimum(Z_vals)
    
    return Ztop, Zbot, Th, La
end

"""
    AdaptiveGrid(X0, Y0, depth, omegaX, omegaY, omegaZ, ax, ay, az, Cr, Nmax)

Calculate the coordinates and respective weights of the point CDMs that form the solution.
Uses the adaptive algorithm explained in reference paper #1.

# Arguments
- `X0`, `Y0`: Horizontal center coordinates
- `depth`: Depth to center
- `omegaX`, `omegaY`, `omegaZ`: Rotation angles (degrees)
- `ax`, `ay`, `az`: Semi-axes
- `Cr`: Grid spacing parameter
- `Nmax`: Maximum number of point CDMs

# Returns
- `Xs`, `Ys`, `Zs`: Coordinates of point CDMs
- `Ws`: Weights of point CDMs
"""
function AdaptiveGrid(X0, Y0, depth, omegaX, omegaY, omegaZ,
                     ax, ay, az, Cr, Nmax)
    
    # Rotation matrix
    Rx = [1  0            0;
          0  cosd(omegaX)  sind(omegaX);
          0 -sind(omegaX)  cosd(omegaX)]
    
    Ry = [cosd(omegaY)  0  -sind(omegaY);
          0             1   0;
          sind(omegaY)  0   cosd(omegaY)]
    
    Rz = [cosd(omegaZ)   sind(omegaZ)  0;
         -sind(omegaZ)   cosd(omegaZ)  0;
          0              0             1]
    
    R = Rz * Ry * Rx
    
    # Determine aC, dT, aV, np
    aC = maximum([ax, ay, az])  # equation 2
    Ztop, Zbot = EllTopBot(depth, omegaX, omegaY, omegaZ, ax, ay, az)
    aV = (Ztop - Zbot) / 2  # equation 4
    dT = -Ztop
    k1 = dT / (aC * Cr - aV)  # equation 6
    
    # Determine aH
    normVec = [0.0, 0.0, 1.0]  # Unit normal vector of a horizontal plane
    normVecR = R' * normVec  # rotated normal vector
    Pt = [0.0, 0.0, 0.0]  # A point on the plane
    PtR = R' * Pt  # A point on the rotated plane
    Pc_tmp, aH = Ell_Plane_Intersect(ax, ay, az, normVecR[1], normVecR[2],
                                     normVecR[3], normVecR' * PtR)
    
    SFa = 0.0
    SFb = 1.0
    Nta = 0
    Ntb = Int(Nmax)
    itr = 1
    
    # Initialize variables that will be set in the loop
    # (needed in case of early break)
    Npar = 0
    Zpar = Float64[]
    Z0 = Float64[]
    Pc = zeros(0, 3)
    ajMax = Float64[]
    ajMin = Float64[]
    ejMax = zeros(0, 3)  # Matrix: each row is a basis vector
    ejMin = zeros(0, 3)  # Matrix: each row is a basis vector
    sjMax = Float64[]
    njMax = Int[]
    sjMin = Float64[]
    njMin = Int[]
    Nsp = Int[]
    Nt = 0
    
    while Nta != Ntb && abs(SFa - SFb) > 1e-3
        
        if itr == 1
            SF = 1.0
        else
            SF = (SFa + SFb) / 2
        end
        
        np = SF * dT ./ (2 * k1 * aV)  # equation 8
        nH = SF * dT ./ (2 * k1 * aH)
        
        # Determine the partitioning planes: Zpar contains the Z coordinates
        Npar0 = floor(Int, log(Zbot / Ztop) / log(1 + 1 / np))  # Initial # of partitions
        
        if Npar0 <= 1
            Zpar = [Ztop; Zbot]
        else
            ZparTmp = Ztop .* (1 + 1/np).^(0:Npar0)
            if ZparTmp[end] > Zbot
                Zpar = [ZparTmp; Zbot]
            else
                Zpar = ZparTmp
            end
        end
        
        Npar = length(Zpar) - 1  # Total number of partitioning planes
        Z0 = (Zpar[1:Npar] .+ Zpar[2:Npar+1]) ./ 2  # Planes containing the pCDMs
        
        # Determine the intersection ellipses
        Ptj = [zeros(2, Npar); (Z0 .+ depth)']  # A point on each plane
        PtjR = R' * Ptj  # A point on the rotated plane
        Pc, ajMax, ajMin, ejMax, ejMin = Ell_Plane_Intersect(ax, ay, az,
            normVecR[1], normVecR[2], normVecR[3], normVecR' * PtjR)
        
        sjMax = -Zpar[1:Npar] ./ nH  # source spacing on the partitioning planes
        njMax = floor.(Int, ajMax ./ sjMax)
        
        sjMin = ajMin ./ ajMax .* sjMax
        njMin = ceil.(Int, ajMin ./ sjMin)
        
        Nsp = zeros(Int, Npar)
        for kk in 1:Npar
            xs, ys = meshgrid(-njMax[kk]:njMax[kk], -njMin[kk]:njMin[kk])
            xs = xs .* sjMax[kk]
            ys = ys .* sjMin[kk]
            Ind = xs.^2 ./ ajMax[kk]^2 .+ ys.^2 ./ ajMin[kk]^2 .< 1
            Nsp[kk] = count(Ind)
        end
        Nt = sum(Nsp)
        
        if itr == 1 && Nt < Nmax
            break
        else
            itr = typemax(Int)  # inf equivalent
        end
        
        if sign(Nt - Nmax) == sign(Nta - Nmax)
            SFa = SF
            Nta = Nt
        else
            SFb = SF
            Ntb = Nt
        end
    end
    
    # Calculate volume fraction of all partitions: Vf
    if Ztop == Zbot
        Vf = [1.0]
    else
        Vp = zeros(Npar + 1)  # Volume of the top ellipsoidal caps
        for kk in 1:(Npar + 1)
            Vp[kk] = Ell_Cap_Vol(depth, omegaX, omegaY, omegaZ, ax, ay, az, Zpar[kk])
        end
        Vp = diff(Vp)  # Volume of the partitions
        Vf = Vp ./ sum(Vp)  # Volume fraction of the partitions
    end
    
    Xs = zeros(Nt)
    Ys = zeros(Nt)
    Zs = zeros(Nt)
    Ws = zeros(Nt)
    Ns0 = 0
    
    for kk in 1:Npar
        
        xs, ys = meshgrid(-njMax[kk]:njMax[kk], -njMin[kk]:njMin[kk])
        xs = xs .* sjMax[kk]
        ys = ys .* sjMin[kk]
        Ind = xs.^2 ./ ajMax[kk]^2 .+ ys.^2 ./ ajMin[kk]^2 .< 1
        xs = xs[Ind]
        ys = ys[Ind]
        
        # Stack basis vectors and coordinates
        # ejMax[kk,:] and ejMin[kk,:] are row vectors, need to make them columns
        basis = hcat(vec(ejMax[kk,:]), vec(ejMin[kk,:]), normVecR)
        coords = [xs'; ys'; zeros(1, length(xs))]
        r = basis * coords
        
        Xr = r[1,:] .+ Pc[kk,1]
        Yr = r[2,:] .+ Pc[kk,2]
        Zr = r[3,:] .+ Pc[kk,3]
        
        r1 = R * [Xr'; Yr'; Zr']
        Xs0 = vec(r1[1,:] .+ X0)
        Ys0 = vec(r1[2,:] .+ Y0)
        Zs0 = vec(r1[3,:] .- depth)
        
        Nsp_local = length(Xs0)
        
        if Npar > 1 && kk == 1
            _, _, _, _, Zs1, Zs2 = Ell_Line_Intersect(X0, Y0, depth,
                omegaX, omegaY, omegaZ, ax, ay, az, 0, 0, 1, Xs0, Ys0, Zs0)
            Wz = vec(maximum(vcat(Zs1, Zs2), dims=1)) .- Zs0[1]
            Ws0 = Wz .+ sjMax[kk] / 2
            sum_Ws0 = sum(Ws0)
            Ws0 = (Vf[kk] / sum_Ws0) .* Ws0
        elseif Npar > 1 && kk == Npar
            _, _, _, _, Zs1, Zs2 = Ell_Line_Intersect(X0, Y0, depth,
                omegaX, omegaY, omegaZ, ax, ay, az, 0, 0, 1, Xs0, Ys0, Zs0)
            Wz = Zs0[1] .- vec(minimum(vcat(Zs1, Zs2), dims=1))
            Ws0 = Wz .+ sjMax[kk] / 2
            sum_Ws0 = sum(Ws0)
            Ws0 = (Vf[kk] / sum_Ws0) .* Ws0
        elseif Npar == 1
            _, _, _, _, Zs1, Zs2 = Ell_Line_Intersect(X0, Y0, depth,
                omegaX, omegaY, omegaZ, ax, ay, az, 0, 0, 1, Xs0, Ys0, Zs0)
            Wz = vec(abs.(Zs1 .- Zs2))
            Ws0 = Wz
            sum_Ws0 = sum(Ws0)
            Ws0 = (Vf[kk] / sum_Ws0) .* Ws0
        else
            Ws0 = fill(Vf[kk] / Nsp_local, Nsp_local)
        end
        
        Xs[Ns0+1:Ns0+Nsp_local] = Xs0
        Ys[Ns0+1:Ns0+Nsp_local] = Ys0
        Zs[Ns0+1:Ns0+Nsp_local] = Zs0
        Ws[Ns0+1:Ns0+Nsp_local] = Ws0
        
        Ns0 = Ns0 + Nsp_local
    end
    
    return Xs, Ys, Zs, Ws
end


"""
    meshgrid(x, y)

Create 2D grid from coordinate vectors (like MATLAB's meshgrid).

# Returns
- `X`, `Y`: 2D coordinate matrices
"""
function meshgrid(x, y)
    X = repeat(reshape(x, 1, :), length(y), 1)
    Y = repeat(y, 1, length(x))
    return X, Y
end

"""
    Ell_Cap_Vol(depth, omegaX, omegaY, omegaZ, ax, ay, az, Z0)

Calculate the ellipsoidal cap volume formed by the intersection of the
horizontal plane z = Z0 and an ellipsoid.

The ellipsoid equation before applying rotations is: x²/ax² + y²/ay² + z²/az² = 1
The plane equation is: z = Z0

# Arguments
- `depth`: Depth to ellipsoid center
- `omegaX`, `omegaY`, `omegaZ`: Rotation angles (degrees)
- `ax`, `ay`, `az`: Semi-axes of ellipsoid
- `Z0`: Z-coordinate of cutting plane

# Returns
- `Vtop`: Volume of ellipsoidal cap above the plane
"""
function Ell_Cap_Vol(depth, omegaX, omegaY, omegaZ, ax, ay, az, Z0)
    
    # Build rotation matrices
    Rx = [1  0            0;
          0  cosd(omegaX)  sind(omegaX);
          0 -sind(omegaX)  cosd(omegaX)]
    
    Ry = [cosd(omegaY)  0  -sind(omegaY);
          0             1   0;
          sind(omegaY)  0   cosd(omegaY)]
    
    Rz = [cosd(omegaZ)   sind(omegaZ)  0;
         -sind(omegaZ)   cosd(omegaZ)  0;
          0              0             1]
    
    R = Rz * Ry * Rx
    
    # Unit normal vector of the plane
    normVec = [0.0, 0.0, 1.0]
    normVecR = R' * normVec  # rotated normal vector
    
    # A point on the plane
    Pt = [0.0, 0.0, Z0 + depth]
    PtR = R' * Pt  # A point on the rotated plane
    
    A = normVecR[1]
    B = normVecR[2]
    C = normVecR[3]
    D = normVecR' * PtR
    
    dn = D / sqrt((A * ax)^2 + (B * ay)^2 + (C * az)^2)
    Vtop = π * ax * ay * az / 3 * (1 - dn)^2 * (2 + dn)
    
    return Vtop
end


"""
    Ell_Plane_Intersect(ax, ay, az, a, b, c, d)

Find the ellipse formed by the intersection of an ellipsoid and a plane.

The ellipsoid equation before rotations: x²/ax² + y²/ay² + z²/az² = 1
The plane equation: ax + by + cz = d

Based on Reference #6: Klein (2012).

# Arguments
- `ax`, `ay`, `az`: Semi-axes of ellipsoid
- `a`, `b`, `c`: Plane normal vector components
- `d`: Distance parameter (can be a vector for multiple planes)

# Returns
- `Pc`: Center points of intersection ellipses (Nx3 matrix if d is vector)
- `aMaj`: Major semi-axes
- `aMin`: Minor semi-axes
- `eMaj`: Major axis unit vectors (3-element vector or 3xN matrix)
- `eMin`: Minor axis unit vectors (3-element vector or 3xN matrix)
"""
function Ell_Plane_Intersect(ax, ay, az, a, b, c, d)
    
    # Find a vector that is normal to [a b c]
    mx, Ind = findmax(abs.([a, b, c]))
    
    if Ind == 1
        r = cross([a, b, c], [0, mx, mx])
    elseif Ind == 2
        r = cross([a, b, c], [mx, 0, mx])
    elseif Ind == 3
        r = cross([a, b, c], [mx, mx, 0])
    end
    
    r = r / sqrt(r' * r)
    normVec = [a, b, c] / sqrt(a^2 + b^2 + c^2)
    s = cross(normVec, r)
    
    # Diagonal matrix
    D1 = diagm([1/ax, 1/ay, 1/az])
    
    w = atan(2 * (D1 * r)' * (D1 * s), 
             ((D1 * r)' * (D1 * r) - (D1 * s)' * (D1 * s))) / 2
    
    e1 = cos(w) * r + sin(w) * s
    e2 = -sin(w) * r + cos(w) * s
    
    # Handle d as vector
    d = d isa AbstractArray ? vec(d) : [d]
    kappa = d ./ sqrt(a^2 + b^2 + c^2)  # distance from the origin
    dn = kappa.^2 .* (a^2 + b^2 + c^2) ./ ((a * ax)^2 + (b * ay)^2 + (c * az)^2)
    
    beta1 = (D1 * e1)' * (D1 * e1)
    beta2 = (D1 * e2)' * (D1 * e2)
    
    a1 = sqrt.((1 .- dn) ./ beta1)
    a2 = sqrt.((1 .- dn) ./ beta2)
    
    if a1[1] > a2[1]
        aMaj = a1
        eMaj = e1
        aMin = a2
        eMin = e2
    else
        aMaj = a2
        eMaj = e2
        aMin = a1
        eMin = e1
    end
    
    Ck = kappa .* (a^2 + b^2 + c^2) ./ ((a * ax)^2 + (b * ay)^2 + (c * az)^2) ./ 
         sqrt(a^2 + b^2 + c^2)
    
    # Build Pc matrix - each row is a center point
    Pc = hcat(Ck .* ax^2 * a, Ck .* ay^2 * b, Ck .* az^2 * c)
    
    # If multiple planes, replicate the basis vectors
    # eMaj and eMin are the same for all planes (only scaling changes)
    if length(d) > 1
        # Create matrices where each row is the same basis vector
        eMaj_mat = repeat(eMaj', length(d), 1)
        eMin_mat = repeat(eMin', length(d), 1)
        return Pc, aMaj, aMin, eMaj_mat, eMin_mat
    else
        return Pc, aMaj, aMin, eMaj, eMin
    end
end


"""
    Ell_Line_Intersect(X0, Y0, depth, omegaX, omegaY, omegaZ, ax, ay, az, a0, b0, c0, x0, y0, z0)

Find the intersection points of a line and an ellipsoid.

Can handle one ellipsoid and multiple lines.
Ellipsoid equation before rotations: x²/ax² + y²/ay² + z²/az² = 1
Line equation: x = k*a0 + x0, y = k*b0 + y0, z = k*c0 + z0

# Arguments
- `X0`, `Y0`, `depth`: Ellipsoid center coordinates
- `omegaX`, `omegaY`, `omegaZ`: Rotation angles (degrees)
- `ax`, `ay`, `az`: Semi-axes
- `a0`, `b0`, `c0`: Line direction components (can be vectors)
- `x0`, `y0`, `z0`: Points on the lines (can be vectors)

# Returns
- `Xs1`, `Xs2`: X-coordinates of intersection points
- `Ys1`, `Ys2`: Y-coordinates of intersection points
- `Zs1`, `Zs2`: Z-coordinates of intersection points
"""
function Ell_Line_Intersect(X0, Y0, depth, omegaX, omegaY, omegaZ, 
                           ax, ay, az, a0, b0, c0, x0, y0, z0)
    
    # Build rotation matrices
    Rx = [1  0            0;
          0  cosd(omegaX)  sind(omegaX);
          0 -sind(omegaX)  cosd(omegaX)]
    
    Ry = [cosd(omegaY)  0  -sind(omegaY);
          0             1   0;
          sind(omegaY)  0   cosd(omegaY)]
    
    Rz = [cosd(omegaZ)   sind(omegaZ)  0;
         -sind(omegaZ)   cosd(omegaZ)  0;
          0              0             1]
    
    R = Rz * Ry * Rx
    
    # Handle scalar direction with vector positions
    if length(a0) == 1 && length(x0) > 1
        a0 = fill(a0, size(x0))
        b0 = fill(b0, size(x0))
        c0 = fill(c0, size(x0))
    end
    
    # Unit vector along the line
    unitVec = [vec_safe(a0)'; vec_safe(b0)'; vec_safe(c0)']
    
    # A point on the line
    P0 = [vec_safe(x0)' .- X0; vec_safe(y0)' .- Y0; vec_safe(z0)' .+ depth]
    
    unitVecR = R' * unitVec  # rotated unit vector
    P0R = R' * P0  # rotated point
    
    x1 = P0R[1, :]
    y1 = P0R[2, :]
    z1 = P0R[3, :]
    
    a = unitVecR[1, :]
    b = unitVecR[2, :]
    c = unitVecR[3, :]
    
    # Intersection of a generic line and a standard ellipsoid
    A = a.^2 .* ay^2 * az^2 .+ b.^2 .* ax^2 * az^2 .+ c.^2 .* ax^2 * ay^2
    B = 2 .* (a .* x1 .* ay^2 * az^2 .+ b .* y1 .* ax^2 * az^2 .+ 
              c .* z1 .* ax^2 * ay^2)
    C = x1.^2 .* ay^2 * az^2 .+ y1.^2 .* ax^2 * az^2 .+ z1.^2 .* ax^2 * ay^2 .- 
        ax^2 * ay^2 * az^2
    
    k1, k2 = SolveQuadraticVec(A, B, C)
    
    xs1 = k1 .* a .+ x1
    xs2 = k2 .* a .+ x1
    ys1 = k1 .* b .+ y1
    ys2 = k2 .* b .+ y1
    zs1 = k1 .* c .+ z1
    zs2 = k2 .* c .+ z1
    
    # Transform back to original coordinate system
    r1 = R * [xs1'; ys1'; zs1']
    r2 = R * [xs2'; ys2'; zs2']
    
    Xs1 = r1[1, :]' .+ X0
    Ys1 = r1[2, :]' .+ Y0
    Zs1 = r1[3, :]' .- depth
    Xs2 = r2[1, :]' .+ X0
    Ys2 = r2[2, :]' .+ Y0
    Zs2 = r2[3, :]' .- depth
    
    return Xs1, Xs2, Ys1, Ys2, Zs1, Zs2
end


"""
    SolveQuadraticVec(A, B, C)

Find the roots of A(k)*x² + B(k)*x + C(k) = 0, where k=1,...,N.

Uses numerically stable formulation to avoid cancellation errors.
A, B, and C must have the same size.

# Arguments
- `A`, `B`, `C`: Coefficient arrays (same size)

# Returns
- `x1`, `x2`: Root arrays (same size as inputs)
"""
function SolveQuadraticVec(A, B, C)
    
    D = B.^2 .- 4 .* A .* C
    
    x1 = zeros(size(A))
    x2 = zeros(size(A))
    
    Ind = B .>= 0
    
    # Numerically stable formulation
    x1[Ind] = (-B[Ind] .- sqrt.(D[Ind])) ./ (2 .* A[Ind])
    x2[Ind] = 2 .* C[Ind] ./ (-B[Ind] .- sqrt.(D[Ind]))
    x1[.!Ind] = 2 .* C[.!Ind] ./ (-B[.!Ind] .+ sqrt.(D[.!Ind]))
    x2[.!Ind] = (-B[.!Ind] .+ sqrt.(D[.!Ind])) ./ (2 .* A[.!Ind])
    
    # In this particular problem D cannot be negative!
    x1[D .< 0] .= 0
    x2[D .< 0] .= 0
    
    return x1, x2
end
