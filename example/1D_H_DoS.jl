using Plots
using LaTeXStrings

include("../src/struct.jl")

# set atom parameters
X1 = [0.0]
X2 = [0.0]
z = 0.0
s1= 1
s2= 1 
L1= √5-1 
L2= 2
v1(x)=0.5*x^2     
v2(x)= 0.5*x^2  
gamma1=2*π/L1
gamma2=2*π/L2
atoms = BiLayer1D(L1=L1, L2=L2, X1=X1, X2=X2, z=z, v1=v1, v2=v2)

# set the model parameters 
EcL =100.0
EcW =20.0
nk = 1
neigs=10
nf= 2 * floor( Int, 2*EcL *L1/π )  
γ=0.01 
h = (2.0*π/L1) / nk
kpoint = range(0.0, 2*π/L1-h, length=nk)
kpts = [kpoint[i] for i=1:length(kpoint)]
model = pwIncommensurate1D_LW(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ);

# Generagte Hamiltonian 
include("../src/1Dincommensurate.jl")

H ,G ,R = hamiltonian(atoms, model);

# test function
g_fermi_beta02mu15(x)= x * (1 + 0.2*exp((x-10)*1))^(-1)
g_fermi_beta1mu10(x)= x * (1 + exp((x-10)*1))^(-1)
g_fermi_beta10mu10(x) = x*(1+exp((x-10)*10))^(-1)
g(x)=g_fermi_beta1mu10(x)

# hamiltonian
H ,G ,R=hamiltonian(atoms, model);
λ = zeros(Float64, size(H,1));
D, Ψ = eigen(Array(H));

# DoS for test function with β=1,μ=10 
d=zeros(Float64,length(D))
for i=1:length(D)
    d[i]=g(D[i])
end
dos=(gamma1*gamma2/EcL)*sum(d)

