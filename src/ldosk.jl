"""Computet the local density of states: g(H(ξ))_{0,0}, g is the gauss function"""


# U0 only include the (0,0) element of eigenvector for all k point
function Ldos(kpts::Vector{Float64}, E::Matrix{Float64}, U0::Array{Float64, 2}; xr = 100, σ = 0.1, h = 0.01)
    λ = real.(E)
    xl = minimum(λ[1, :]) - 0.1
    #xr = 100 #maximum(λ[end, :]) + 0.1
    x = range(xl, stop=xr, length=floor(Int, (xr - xl) / h))
    LDos_k = zeros(Float64, length(x), length(kpts))
    E_x = collect(x)
    for k = 1 : length(kpts)
        y = zeros(Float64, length(x))
        for i = 1 : length(x), j = 1 : size(λ, 1)
            y[i] = y[i] + exp(-((x[i] - λ[j, k]) / σ)^2 / 2) / (σ * √(2 * π)) * norm(U0[j, k])^2 #/Lz
        end
        LDos_k[:, k] = y
    end
    return E_x, LDos_k 
end