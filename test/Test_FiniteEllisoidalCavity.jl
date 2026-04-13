using Test

include("../src/FiniteEllipsoidalCavity.jl")

# make mesh
x = 0:100:20e3
y = 0:200:20e3
X, Y = meshgrid(x, y)
xl = extrema(x) ./ 1e3; yl = extrema(y) ./ 1e3

# cavity properties
X0       = 0
Y0       = 0
depth    = 10250
DepthRef = "C"
omegaX   = 0
omegaY   = 0
omegaZ   = 0
ax       = 550
ay       = 550
az       = 3750

# pressure
p        = 18.8e6

# crust elasticity
nu       = 0.25
mu       = 10e9
lambda   = (2*mu*nu)/(1-2*nu)

# numerics
Nmax     = 5e3
Cr       = 14

# compute
ue,un,uv,dV,DV,Ns = fECM(X,Y,X0,Y0,depth,omegaX,omegaY,omegaZ,ax,ay,az,p,mu,lambda,DepthRef;Nmax=Nmax,Cr=Cr)

# compare some results to original Matlab version
@test all(isapprox.(extrema(uv), (8.044701021001343e-04, 0.005989663901181), rtol=1e-4))
@test uv[1,1]          ≈ 0.004201484984743  rtol=1e-4
@test uv[9,100]        ≈ 0.004833273685574  rtol=1e-4
@test uv[87,46]        ≈ 0.002210474017935  rtol=1e-4
@test maximum(ue)      ≈ 0.005037030511448  rtol=1e-4
@test maximum(un)      ≈ 0.005037030511448  rtol=1e-4

# change some parameters
depth    = 5000
omegaY   = 25
omegaZ   = -65
ax       = 850
az       = 1500

# compute again
ue,un,uv,dV,DV,Ns = fECM(X,Y,X0,Y0,depth,omegaX,omegaY,omegaZ,ax,ay,az,p,mu,lambda,DepthRef;Nmax=Nmax,Cr=Cr)

# compare some results to original Matlab version
@test all(isapprox.(extrema(uv), (3.405925303200454e-04, 0.022797797325954), rtol=1e-4))
@test uv[1,1]          ≈ 0.020284896020590  rtol=1e-4
@test uv[9,100]        ≈ 0.004740361583832  rtol=1e-4
@test uv[87,46]        ≈ 0.001025808009311  rtol=1e-4
@test maximum(ue)      ≈ 0.017285000082854  rtol=1e-4
@test maximum(un)      ≈ 0.011987093395042  rtol=1e-4

