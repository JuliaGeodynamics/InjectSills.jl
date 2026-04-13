"""
    Mogi_McTigue(dP, nu, G, r, d, x)

Computes the surface displacement due to a pressurized spherical cavity according to Mogi, 1958 and McTigue, 1987. Units of output will be consistent with units of input. 

# Arguments
- `dP`: Pressure change in cavity
- `nu`: Poisson's ratio of the host rock
- `G`:  Shear modulus of the host rock
- `r`:  Radius of the cavity
- `d`:  Depth of the center of the Sphere
- `x`:  Horizontal distance of observation point from center of cavity. Can be a vector or matrix

# Outputs
- `Mogi_Uz`:    Vertical displacement according to Mogi
- `Mogi_Ur`:    Horizontal displacement according to Mogi
- `McTigue_Uz`: Vertical displacement according to McTigue
- `McTigue_Ur`: Horizontal displacement according to McTigue
"""
function Mogi_McTigue(dP, nu, G, r, d, x)
    R = sqrt.(x.^2.0 .+ d.^2.0)
    
    Mogi_Uz = r^3 .* dP .* (1.0-nu) .* d ./ (G .* R.^3)
    Mogi_Ur = r^3 .* dP .* (1.0-nu) .* x ./ (G .* R.^3)
    
    term1   = (r ./ d) .^ 3.0
    term2   = (1.0+nu) / (2.0 * (-7.0 + 5.0 * nu))
    term3   = (15.0 * d.^2.0 * (-2.0 + nu)) ./ (4.0 .* R.^2.0 * (-7.0 + 5.0 * nu))
    
    McTigue_corr = 1.0 .+ term1 .* (term2 .+ term3)
    
    McTigue_Uz = Mogi_Uz .* McTigue_corr
    McTigue_Ur = Mogi_Ur .* McTigue_corr

    return Mogi_Uz, Mogi_Ur, McTigue_Uz, McTigue_Ur
end