using StaticArrays
using SparseArrays
using LinearAlgebra
using Arpack
using FFTW

function hamiltonian(atoms::BiLayer1D, model::pwIncommensurate1D_LW)
    L1 = atoms.L1
    L2 = atoms.L2
    X1 = atoms.X1
    X2 = atoms.X2
    v1 = atoms.v1
    v2 = atoms.v2
    EcL = model.EcL
    EcW = model.EcW
    kpts = model.kpts
    γ = model.γ

    m_max1 = floor(Int, max(EcL, EcW) * L1 / (2.0 * pi))
    m_max2 = floor(Int, max(EcL, EcW) * L2 / (2.0 * pi))

    # calculate the potential and perform fftw
    n_fftw = model.n_fftw
    grids1 = range(0, L1 - L1 / n_fftw, length = n_fftw)
    vreal1 = zeros(Float64, length(grids1))
    for j = 1:length(X1)
        for i = 1:n_fftw
            vreal1[i] += v1(grids1[i] - X1[j])
        end
    end
    vfft1 = fft(vreal1) / n_fftw
    grids2 = range(0, L2 - L2 / n_fftw, length = n_fftw)
    vreal2 = zeros(Float64, length(grids2))
    for j = 1:length(X2)
        for i = 1:n_fftw
            vreal2[i] += v2(grids2[i] - X2[j])
        end
    end
    vfft2 = fft(vreal2) / n_fftw
    # count the degrees of freedom and label the basis functions
    npw = 0
    G = Int[]
    for j1 = -m_max1:m_max1, j2 = -m_max2:m_max2
        G1 = 2.0 * pi * j1 / L1
        G2 = 2.0 * pi * j2 / L2
        #if norm([G1+G2]) < EcW
        #if norm([G1-G2])<EcL
        if abs(G1 + G2) < EcW && abs(G1 - G2) < EcL
            #if norm(G1,G2)^2 < 2.0 * Ec    #if norm([G1,G2])^2 < 2.0 * Ec
            npw = npw + 1    # all coupled points statisfied Ec
            append!(G, [j1, j2])
        end
    end
    G = reshape(G, 2, Int(length(G) / 2))'


    println("EcL = ", EcL, "EcW = ", EcW, "; DOF = ", npw)

    # assemble the matrix
    #nk = length(kpts)
    val = Complex{Float64}[]
    indi = Int64[]
    indj = Int64[]
    R = zeros(Float64, npw)
    diffG = zeros(Int64, 2)
    for n1 = 1:npw, n2 = 1:npw
        G11 = 2.0 * pi * G[n1, 1] / L1 # n1th Gm
        G12 = 2.0 * pi * G[n1, 2] / L2 # n1th Gn
        G21 = 2.0 * pi * G[n2, 1] / L1
        G22 = 2.0 * pi * G[n2, 2] / L2
        #diffG = G[n1] - G[n2]
        @views diffG = G[n1, :] - G[n2, :]
        cp = (
            n1 == n2 ?
            0.5 * (G11^2 + G12^2) + G11 * G12
            : (G[n1, 2] == G[n2, 2] ? exp(-γ * (2.0 * pi * diffG[1] / L1)^2) : 0.0)
            +
            (G[n1, 1] == G[n2, 1] ? exp(-γ * (2.0 * pi * diffG[2] / L2)^2) : 0.0)
        )

        if cp != 0.0
            push!(val, cp)
            push!(indi, n1)
            push!(indj, n2)
        end

        R[n1] = G11 + G12
    end
    H = sparse(indi, indj, val)

    return H, G, R 
end
