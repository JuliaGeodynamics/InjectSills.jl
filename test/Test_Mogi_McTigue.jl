using Test

include("../src/Mogi_McTigue.jl")

# make coordinates
x        = collect(0:100:10e3)

# cavity properties
d        = 5000
r        = 1500
dP       = 10e6

# crust properties
G        = 10e9
nu       = 0.25

# compute
Mogi_Uz, Mogi_Ur, McTigue_Uz, McTigue_Ur = Mogi_McTigue(dP, nu, G, r, d, x)

# check results
@test Mogi_Uz[1]     ≈ 0.10125
@test McTigue_Uz[1]  ≈ 0.10407289402173914
@test Mogi_Uz[36]    ≈ 0.05566928318963278
@test McTigue_Uz[36] ≈ 0.05665722209549374
@test Mogi_Ur[1]     ≈ 0.0
@test McTigue_Ur[1]  ≈ 0.0
@test Mogi_Ur[36]    ≈ 0.038968498232742954
@test McTigue_Ur[36] ≈ 0.039660055466845624