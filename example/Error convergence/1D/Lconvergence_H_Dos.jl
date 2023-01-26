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

# Generate EcL 
Lstep = 0.05
lengthL = 25
vecEcL = zeros(Float64,lengthL);
initialL = 2
for i = 1 : lengthL  
    vecEcL[i] = 10^(initialL + Lstep*(i-1))
end

# test function
g_fermi_beta02mu15(x) = x*(1+exp((x-15)*0.2))^(-1)
g(x) = g_fermi_beta02mu15(x)

dos_fermi_beta02mu15=zeros(Float64,lengthL);
for i = 1 : lengthL
    EcL = vecEcL[i]
    EcW = 100.0
    neigs = 10
    nf = 2 * floor( Int, 2*EcL *L1/π )  
    γ = 0.05 
    kpts = [0.]
    model = pwIncommensurate1D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ)
    # hamiltonian
    H ,G ,R = hamiltonian(atoms, model);
    λ = zeros(Float64, size(H,1));
    D, Ψ = eigen(Array(H));
    d_fermi_beta02mu15_EcW = zeros(Float64,length(D))
    for n = 1 : length(D)
        d_fermi_beta02mu15_EcW[n] = g(D[n]) 
    end
    dos_fermi_beta02mu15[i]=(gamma1*gamma2/EcL)*sum(d_fermi_beta02mu15_EcW)
end

# Exact solution EcL = 5000, EcW = 100
EcL = 5000
model = pwIncommensurate1D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ)
H ,G ,R = hamiltonian(atoms, model);
λ = zeros(Float64, size(H,1));
D, Ψ = eigen(Array(H));
d_fermi_beta02mu15_EcW = zeros(Float64,length(D))
for n = 1 : length(D)
    d_fermi_beta02mu15_EcW[n] = g(D[n]) 
end
dos_fermi_beta02mu15_exact = (gamma1*gamma2/EcL)*sum(d_fermi_beta02mu15_EcW)

# plot the error
Lerror_fermi_beta02mu15 = zeros(Float64,lengthL)
for i = 1 : lengthL
    Lerror_fermi_beta02mu15[i] = abs(dos_fermi_beta02mu15_exact - dos_fermi_beta02mu15[i])
end   
plot(vecEcL, Lerror_fermi_beta02mu15, xaxis=:log, yaxis=:log, label="Error", xlabel="L", ylabel= "Error", ls=:solid, legend=:bottomleft, shape=:circle, markersize=6, lw=3, tickfontsize=15, legendfontsize=15, guidefontsize=15)
plot!(vecEcL, 20*vecEcL.^(-1), xaxis=:log, yaxis=:log, label= L"\mathcal{O}(L^{-1})", ls=:dash, lw=3, c=:black)