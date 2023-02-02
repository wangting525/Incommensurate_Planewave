"""`dos`: the density of states of the given set of eigenvalues
    λ : eigenvalues
    σ : width for smearing
    h : meshsize for the density of states
"""
function dos(λ::Vector{Float64}, σ::Float64, h::Float64)
    ϵ = λ[:]
    xl = minimum(ϵ[:]) - 2 * σ
    xr = maximum(ϵ[:]) + 2 * σ
    x = range(xl, stop=xr, length=floor(Int, (xr - xl) / h))
    y = zeros(length(x), 1)
    for i = 1:length(x), j = 1:length(ϵ)
        y[i] = y[i] + exp(-((x[i] - ϵ[j]) / σ)^2 / 2) / (σ * √(2 * π))
    end
    #plot(x, y)
    return x, y
end
