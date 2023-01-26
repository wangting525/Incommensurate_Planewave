"""W convergence for different γ,"""
using Plots
using LaTeXStrings

include("../../../src/struct.jl")
include("../../../src/1Dincommensurate.jl")

X1 = [0.0]
X2 = [0.0]
z = 0.0
s1= 1
s2= 1 
L1= √5-1 
L2= 2
gamma1=2*π/L1
gamma2=2*π/L2
atoms = BiLayer1D(L1=L1, L2=L2, X1=X1, X2=X2, z=z)

g_fermi_beta1mu10(x) = x*(1+exp((x-10)*1))^(-1)
g(x) = g_fermi_beta1mu10(x)

# Generate EcW 
Wstep = 5
lengthW = 7
vecEcW = zeros(Int64,lengthW);
intialW = 5;
for i=1:lengthW
    vecEcW[i]=intialW+Wstep*(i-1)
end

"""γ = 0.001"""
dos_fermi_beta1mu10_gamma0001 = zeros(Float64,lengthW);
for i=1:lengthW
    EcL =2000
    EcW =vecEcW[i]
    nk = 1
    neigs=10
    nf= 2 * floor( Int, 2*EcL *L1/π )  
    γ=0.001 
    kpts = [0.]

    model = pwIncommensurate1D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ)
    # hamiltonian
    H ,G ,R=hamiltonian(atoms, model);
    λ = zeros(Float64, size(H,1));
    D, Ψ = eigen(Array(H));
    d_fermi_beta1mu10_EcW = zeros(Float64,length(D))
    for n=1:length(D)
        d_fermi_beta1mu10_EcW[n]=g(D[n])
    end
    dos_fermi_beta1mu10_gamma0001[i]=(gamma1*gamma2/EcL)*sum(d_fermi_beta1mu10_EcW)
end

# exact solution 
EcW =100
γ=0.001  # initial state coefficients  ϕ0[i,j]=exp(-γ*(norm(R[i,j]))^2), and R[i,j]=G1[i]+G2[j]
kpts = [0.]
model = pwIncommensurate1D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ)
# hamiltonian
H ,G ,R=hamiltonian(atoms, model);
λ = zeros(Float64, size(H,1));
D, Ψ = eigen(H);
d_fermi_beta1mu10_gamma0001_exact = zeros(Float64,length(D))
for i=1:length(D)
	d_fermi_beta1mu10_gamma0001_exact[i]=g(D[i])
end
dos_fermi_beta1mu10_gamma0001_exact=(gamma1*gamma2/EcL)*sum(d_fermi_beta1mu10_gamma0001_exact)
# Error for gamma = 0.001
Werror


"""γ = 0.01"""
dos_fermi_beta1mu10_gamma001 = zeros(Float64,lengthW);
for i=1:lengthW
    EcL =2000
    EcW =vecEcW[i]
    nk = 1
    neigs=10
    nf= 2 * floor( Int, 2*EcL *L1/π )  
    γ=0.01 
    kpts = [0.]
    model = pwIncommensurate1D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ)
    # hamiltonian
    H ,G ,R=hamiltonian(atoms, model);
    λ = zeros(Float64, size(H,1));
    D, Ψ = eigen(Array(H));
    d_fermi_beta1mu10_EcW = zeros(Float64,length(D))
    for n=1:length(D)
        d_fermi_beta1mu10_EcW[n]=g(D[n])
    end
    dos_fermi_beta1mu10_gamma001[i]=(gamma1*gamma2/EcL)*sum(d_fermi_beta1mu10_EcW)
end

# exact solution 
EcW =100
γ=0.01  # initial state coefficients  ϕ0[i,j]=exp(-γ*(norm(R[i,j]))^2), and R[i,j]=G1[i]+G2[j]
kpts = [0.]
model = pwIncommensurate1D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ)
# hamiltonian
H ,G ,R=hamiltonian(atoms, model);
λ = zeros(Float64, size(H,1));
D, Ψ = eigen(H);
d_fermi_beta1mu10_gamma001_exact = zeros(Float64,length(D))
for i=1:length(D)
	d_fermi_beta1mu10_gamma001_exact[i]=g(D[i])
end
dos_fermi_beta1mu10_gamma001_exact=(gamma1*gamma2/EcL)*sum(d_fermi_beta1mu10_gamma001_exact)


"""γ = 0.1"""
dos_fermi_beta1mu10_gamma01 = zeros(Float64,lengthW);
for i=1:lengthW
    EcL =2000
    EcW =vecEcW[i]
    nk = 1
    neigs=10
    nf= 2 * floor( Int, 2*EcL *L1/π )  
    γ=0.1 
    kpts = [0.]
    model = pwIncommensurate1D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ)
    # hamiltonian
    H ,G ,R=hamiltonian(atoms, model);
    λ = zeros(Float64, size(H,1));
    D, Ψ = eigen(Array(H));
    d_fermi_beta1mu10_EcW = zeros(Float64,length(D))
    for n=1:length(D)
        d_fermi_beta1mu10_EcW[n]=g(D[n])
    end
    dos_fermi_beta1mu10_gamma01[i]=(gamma1*gamma2/EcL)*sum(d_fermi_beta1mu10_EcW)
end

# exact solution 
EcW =100
γ=0.1  # initial state coefficients  ϕ0[i,j]=exp(-γ*(norm(R[i,j]))^2), and R[i,j]=G1[i]+G2[j]
kpts = [0.]
model = pwIncommensurate1D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ)
# hamiltonian
H ,G ,R=hamiltonian(atoms, model);
λ = zeros(Float64, size(H,1));
D, Ψ = eigen(H);
d_fermi_beta1mu10_gamma01_exact = zeros(Float64,length(D))
for i=1:length(D)
	d_fermi_beta1mu10_gamma01_exact[i]=g(D[i])
end
dos_fermi_beta1mu10_gamma01_exact=(gamma1*gamma2/EcL)*sum(d_fermi_beta1mu10_gamma01_exact)

# Error for different γ
Werror_Fermi_beta1_mu_10_gamma0001 = zeros(Float64,lengthW);
Werror_Fermi_beta1_mu_10_gamma001 = zeros(Float64,lengthW);
Werror_Fermi_beta1_mu_10_gamma01 = zeros(Float64,lengthW);
for i = 1 : lengthW
	Werror_Fermi_beta1_mu_10_gamma0001[i]=abs(dos_fermi_beta1mu10_gamma0001_exact - dos_fermi_beta1mu10_gamma0001[i])
	Werror_Fermi_beta1_mu_10_gamma001[i]=abs(dos_fermi_beta1mu10_gamma001_exact - dos_fermi_beta1mu10_gamma001[i])
	Werror_Fermi_beta1_mu_10_gamma01[i]=abs(dos_fermi_beta1mu10_gamma01_exact - dos_fermi_beta1mu10_gamma01[i])
end

Werror = [Werror_Fermi_beta1_mu_10_gamma0001,Werror_Fermi_beta1_mu_10_gamma001,Werror_Fermi_beta1_mu_10_gamma01]
plot(vecEcW, Werror, yaxis=:log, label=[L"\gamma=0.001" L"\gamma=0.01" L"\gamma=0.1"],  
     xlabel="W", ylabel= "Error", legend= :bottomleft,
 shape=[:circle :utriangle :rect], markersize=6, lw=3, tickfontsize=15, legendfontsize=15, guidefontsize=15)
