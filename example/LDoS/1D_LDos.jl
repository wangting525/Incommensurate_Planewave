include("../../src/struct.jl")
include("../../src/1Dincommensurate_shiftedH.jl")
include("../../src/ldosk.jl")

using Plots
using LaTeXStrings
using JLD
using Measures

# set atom parameters
X1 = [0.0]
X2 = [0.0]
z = 0.0
s1= 1 
s2= 1 
L1= √5-1 
L2= 2
atoms = BiLayer1D(L1=L1, L2=L2, X1=X1, X2=X2, z=z)

# set the model parameters 
EcL = 20.0
EcW = 15.0
h = 0.01 
nk = floor(Int, 2 * EcW / h)
neigs = 10
nf= 2 * floor(Int, 2*EcL *L1/π )  
γ = 0.01  
kpoint = range(-EcW, EcW-h, length=nk)
kpts = [kpoint[i] for i=1:length(kpoint)]
model = pwIncommensurate1D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ);

# hamiltonian
H, G, R = hamiltonian(atoms, model);

zero_index = findfirst(x -> x[1] == 0 && x[2] == 0, G)

E = zeros(Float64, size(H,2), nk);
U = zeros(Float64, size(H,2), size(H,2), nk);
for k = 1 : nk
    E[:, k], ψ = eigen(H[k,:,:]); 
	Q, R = qr(ψ)
    U[:,:,k] = Matrix(Q)
end
 
U_0=zeros(Float64, size(H,2), nk);

for k =1 : nk
	U_0[:, k] = U[zero_index, :, k]
end   

# xr is the upper bound of energy, σ is the width of Gauss, h is the meshsize of the energy.
E_x, ldosk = Ldos(kpts, E, U_0; xr= 10, σ = 0.5, h = 0.01)

gr(fmt=:png, dpi=300)
heatmap(kpts, E_x, ldosk, xlabel = L"\xi", xlim = (-8, 8), ylim = (E_x[1], 10), ylabel = "Energy",tickfontsize=15, legendfontsize=15, guidefontsize=15, margin = 2.5mm)