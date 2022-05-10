# module incommensurate2d (for layer)
using StaticArrays
using SparseArrays
using LinearAlgebra
using Arpack
using FFTW

"""`Hamiltonian`:
Hamiltonian for 1D/2D system, using plane wave discretizations for
(-1/2 Δ + V)u = λu
with V = ∑_{k≠0} Z/|k²| ⋅ e^{ikx}/|Γ|^{1/2} satisfying -(Δ+m)V = 4πZδ
basis function: ψ_k = e^{ikx}/|Γ|^{1/2}
"""

function hamiltonian2d(atoms::TwoLayerIn2D, model::pwIncommensurate2D)
    R1 = atoms.R1
    R2 = atoms.R2
    B1 = atoms.B1
    B2 = atoms.B2
    v1 = atoms.v1
    v2 = atoms.v2
    n_fftwx = model.n_fftwx
    n_fftwy = model.n_fftwy

    #Ec = model.Ecut
    kpts = model.kpts
    EcL = model.EcL
    EcW = model.EcW
    γ = model.γ

    gridsx = range(0, L - L / n_fftwx, length = n_fftwx)
    gridsy = range(0, L - L / n_fftwy, length = n_fftwy)
    vreal1 = zeros(Float64, length(gridsx), length(gridsy))
    vreal2 = zeros(Float64, length(gridsx), length(gridsy))
    for i = 1:length(gridsx)
        for j = 1:length(gridsy)
            #lattice1x=A[1,1]*gridsx[i]+A[1,2]*gridsy[j]
            vreal1[i, j] = v1(gridsx[i], gridsy[j])
        end
    end
    for i = 1:length(gridsx)
        for j = 1:length(gridsy)
            #lattice1x=A[1,1]*gridsx[i]+A[1,2]*gridsy[j]
            vreal2[i, j] = v2(gridsx[i], gridsy[j])
        end
    end

    vfft1 = fft(vreal1) / (n_fftwx * n_fftwy)
    vfft2 = fft(vreal2) / (n_fftwx * n_fftwy)

    Gmax11 = floor(Int, max(EcL, EcW) / norm(B1[:, 1]))
    Gmax12 = floor(Int, max(EcL, EcW) / norm(B1[:, 2]))
    Gmax21 = floor(Int, max(EcL, EcW) / norm(B2[:, 1]))
    Gmax22 = floor(Int, max(EcL, EcW) / norm(B2[:, 2]))

    # area of the primative unit cell        a*b*sinθ
    S1 = sqrt(norm(R1[:, 1])^2 * norm(R1[:, 2])^2 - dot(R1[:, 1], R1[:, 2])^2)
    S2 = sqrt(norm(R2[:, 1])^2 * norm(R2[:, 2])^2 - dot(R2[:, 1], R2[:, 2])^2)

    # area of the reciprocal unit cell
    RS1 = sqrt(norm(B1[:, 1])^2 * norm(B1[:, 2])^2 - dot(B1[:, 1], B1[:, 2])^2)
    RS2 = sqrt(norm(B2[:, 1])^2 * norm(B2[:, 2])^2 - dot(B2[:, 1], B2[:, 2])^2)

    # count the degrees of freedom and label the basis functions
    # compute npw

    G = Int[]
    Gmn = Float64[]
    npw = 0
    G1m = zeros(2)
    G2n = zeros(2)
    for j11 = -Gmax11:Gmax11, j12 = -Gmax12:Gmax12
        for j21 = -Gmax21:Gmax21, j22 = -Gmax22:Gmax22
            @views B11 = B1[:, 1]
            @views B12 = B1[:, 2]
            @views B21 = B2[:, 1]
            @views B22 = B2[:, 2]
            @. G1m = j11 * B11 + j12 * B12    #  G1m=[G1_m1;G1_m2]
            @. G2n = j21 * B21 + j22 * B22
            if norm(G1m + G2n) < EcW && norm(G1m - G2n) < EcL
                npw = npw + 1
                append!(G, [j11, j12, j21, j22])
                append!(Gmn, [G1m[1], G1m[2], G2n[1], G2n[2]])
            end
        end
    end
    G = reshape(G, 4, Int(length(G) / 4))'
    Gmn = reshape(Gmn, 4, Int(length(Gmn) / 4))'


    println(" EcutL = ", EcL, "; EcutW = ", EcW, "; DOF = ", npw)


    val = Complex{Float64}[]
    indi = Int64[]
    indj = Int64[]

    R = zeros(Float64, 2, npw)
    diffG = zeros(Int64, 4)

    BG1 = zeros(2)
    BG2 = zeros(2)

    for n1 = 1:npw, n2 = 1:npw

        @views g1 = G[n1, :]
        @views g2 = G[n2, :]
        @. diffG = g1 - g2

        @views gmn = Gmn[n1, :]
        @views d1 = diffG[1:2]
        @views d2 = diffG[3:4]

        cp = (
            n1 == n2 ?
            0.5 * (norm(gmn)^2 + 2 * (gmn[1] * gmn[3] + gmn[2] * gmn[4])) + 1 + 1 :
            #( G[3:4,G1] == G[3:4,G2] ?
            (g1[3] == g2[3] && g1[4] == g2[4] ?
             exp(-γ * norm(mul!(BG1, B1, d1))^2) : 0.)
            #vfft1[ mod((G[1,G1]-G[1,G2]),nfx)+1, mod((G[2,G1]-G[2,G2]),nfy)+1 ] :  0 )
            # * exp( im * dot(G[1:2,G2]-G[1:2,G1], Ratom) )
            #+ ( G[1:2,G1] == G[1:2,G2] ?
            +
            (g1[1] == g2[1] && g1[2] == g2[2] ?
             exp(-γ * norm(mul!(BG2, B2, d2))^2) : 0.)
            #vfft2[ mod((G[3,G1]-G[3,G2]),nfx)+1, mod((G[4,G1]-G[4,G2]),nfy)+1 ] : 0 )
            # * exp( im * dot(G[3:4,G2]-G[3:4,G1], Ratom) )
        )

        if cp != 0.0
            push!(val, cp)
            push!(indi, n1)
            push!(indj, n2)
        end

        R[1, n1] = gmn[1] + gmn[3]
        R[2, n1] = gmn[2] + gmn[4]

    end
    H = sparse(indi, indj, val)

    return H, G, Gmn, R, Gmax11, Gmax12, Gmax21, Gmax22, S1, S2, RS1, RS2
end
