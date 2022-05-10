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
    γ=model.γ
    m_max1 = floor( Int, max(EcL,EcW) * L1 / (2.0 * pi) )
    #m_max1 = floor( Int, EcL * L1 / (2.0 * pi) )
    m_max2 = floor( Int, max(EcL,EcW) * L2 / (2.0 * pi) )
    #m_max2 = floor( Int, EcL * L2 / (2.0 * pi) )
    # calculate the potential and perform fftw
    n_fftw = model.n_fftw
    grids1 = range(0, L1-L1/n_fftw, length = n_fftw)
    vreal1 = zeros(Float64, size(grids1))
    for j = 1 : length(X1)
        vreal1 .+= [ v1(grids1[i]-X1[j]) for i=1:n_fftw ]
    end
    vfft1 = fft(vreal1) / n_fftw
    grids2 = range(0, L2-L2/n_fftw, length = n_fftw)
    vreal2 = zeros(Float64, size(grids2))
    for j = 1 : length(X2)
        vreal2 .+= [ v2(grids2[i]-X2[j]) for i=1:n_fftw ]
    end
    vfft2 = fft(vreal2) / n_fftw
    # count the degrees of freedom and label the basis functions
    npw = 0
    G = []
    #Gmn =[]
    for j1 = -m_max1 : m_max1, j2 = -m_max2 : m_max2
        G1 = 2.0 * pi * j1/L1;
        G2 = 2.0 * pi * j2/L2;
        #if norm([G1+G2]) < EcW
           #if norm([G1-G2])<EcL
        if  abs(G1+G2) < EcW
           if abs(G1-G2) < EcL
           #if norm(G1,G2)^2 < 2.0 * Ec    #if norm([G1,G2])^2 < 2.0 * Ec
              npw = npw + 1    # all coupled points statisfied Ec
              push!(G, [j1 j2]')
              #push!(Gmn,[G1 G2]')
           end
        end
    end

    println("EcL = ", EcL, "EcW = ", EcW, "; DOF = ", npw)
    # assemble the matrix
    nk = length(kpts)
    H = zeros(Complex{Float64}, nk, npw, npw)
    R = zeros(Float64,npw,1)
    #NR = zeros(Float64,nk,npw,1)
    for k = 1 : nk
        println("assembling H[k] with kpt = ", kpts[k])
        for n1 = 1 : npw, n2 = 1 : npw
            G11 = 2.0 * pi * G[n1][1]/L1; # n1th Gm
            G12 = 2.0 * pi * G[n1][2]/L2; # n1th Gn
            G21 = 2.0 * pi * G[n2][1]/L1;
            G22 = 2.0 * pi * G[n2][2]/L2;
            diffG=G[n1][:]-G[n2][:]

            H[k, n1, n2] = ( n1==n2 ?
                0.5 * norm([G11+kpts[k] G12])^2 + (G11+kpts[k]) * G12 + 1 + 1
                : ( G[n1][2]==G[n2][2] ? exp(-γ*norm(2.0 * pi * diffG[1]/L1)^2) : 0 )
                + ( G[n1][1]==G[n2][1] ? exp(-γ*norm(2.0 * pi * diffG[2]/L2)^2) : 0 )
                )
            R[n1, 1]=G11+G12
        end
    end

    return H, G, R 
end
