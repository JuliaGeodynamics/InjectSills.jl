# Numerical checks of displacement fields: quadrature over lines, and stresses by finite
# differences.
using LinearAlgebra

# Gauss–Legendre nodes and weights on [-1, 1] (Golub–Welsch).
function gausslegendre(n)
    β = [k / sqrt(4k^2 - 1) for k in 1:n-1]
    E = eigen(SymTridiagonal(zeros(n), β))
    return E.values, 2 .* E.vectors[1, :] .^ 2
end
const GL_X, GL_W = gausslegendre(40)

# ∫₀^∞ f(s) ds: Gauss–Legendre panels on [0, L] with edges `breaks` refined geometrically
# towards them, and the tail s = L/u, u ∈ (0, 1].
function halfline_integral(f, L; breaks=(), panels=200)
    edges = sort(unique([0.0; L; collect(breaks)]))
    nodes = Float64[]
    for k in 1:length(edges)-1
        a, b = edges[k], edges[k+1]
        # geometric grading towards both ends of [a, b]
        g = [0.5 - 0.5 * cospi(j / panels) for j in 0:panels]
        append!(nodes, a .+ (b - a) .* g)
    end
    nodes = unique(nodes)
    tot = 0.0
    for k in 1:length(nodes)-1
        a, b = nodes[k], nodes[k+1]
        for (x, w) in zip(GL_X, GL_W)
            tot += w * (b - a) / 2 * f((a + b) / 2 + (b - a) / 2 * x)
        end
    end
    for (x, w) in zip(GL_X, GL_W)
        u = (x + 1) / 2
        tot += w / 2 * f(L / u) * L / u^2
    end
    return tot
end

# ∫ f dx over the real line
line_integral(f, L; breaks=()) =
    halfline_integral(f, L; breaks) + halfline_integral(x -> f(-x), L; breaks)

# Displacement gradient J[i, j] = ∂u_i/∂x_j of u(x::Vector) by central differences.
function displacement_gradient(u, x, h)
    n = length(x)
    J = zeros(n, n)
    for j in 1:n
        e = zeros(n); e[j] = h
        J[:, j] = (u(x + e) - u(x - e)) / 2h
    end
    return J
end

# Isotropic stress σ/G from the displacement gradient, with λ/G = 2ν/(1-2ν).
function stress(u, x, h, ν)
    J = displacement_gradient(u, x, h)
    ε = (J + J') / 2
    return 2ν / (1 - 2ν) * tr(ε) * I + 2ε
end
