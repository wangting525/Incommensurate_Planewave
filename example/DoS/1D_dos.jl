using Plots
using LaTeXStrings

include("../../src/struct.jl")
include("../../src/1Dincommensurate.jl")
include("../../src/dos.jl")

# set atom parameters
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

# set the model parameters 
EcL =1000.0
EcW =50.0
nk = 1
neigs=10
nf= 2 * floor( Int, 2*EcL *L1/π )  
γ = 0.01      
kpts = [0.]
model = pwIncommensurate1D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ);

# compute Dos and plot Dos
H, G, R=hamiltonian(atoms, model);
D, Ψ = eigen(Array(H));
λ = real.(D)
σ = 0.1   # width for smearing
h = 0.05 # meshsize for the density of states
x, y = dos(λ, σ, h)
# plot the first 100 points
plot(x[1:100], y[1:100], label = "EcL = 1000, EcW = 50", ylabel="DoS", lw = 3)


# Plot Dos convergence with respect to L
vecEcL = [500, 600, 800, 1000] 
EcW =50.0

X_L = []
DoS_L = []
for i = 1 : length(vecEcL)
	EcL = vecEcL[i]
	model = pwIncommensurate1D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ);
	H, G, R=hamiltonian(atoms, model);
	D, Ψ = eigen(Array(H));
	λ = real.(D)
	σ = 0.1   # width for smearing
	h = 0.05 # meshsize for the density of states
    x, y = dos(λ, σ, h)
    push!(X_L, x[1:100])
	push!(DoS_L, y[1:100])
end
plot(X_L, DoS_L, label = ["EcL = 500, EcW = 50" "EcL = 600, EcW = 50" "EcL = 800, EcW = 50" "EcL = 1000, EcW = 50"], ylabel="DoS", lw = 3)


# plot Dos convergence with respect to W
EcL = 1000 
vecEcW = [5, 10, 50]; 

X_W = []
DoS_W = []
for i = 1 : length(vecEcW)
	EcW = vecEcW[i]
	model = pwIncommensurate1D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ);
	H, G, R=hamiltonian(atoms, model);
	D, Ψ = eigen(Array(H));
	λ = real.(D)
	σ = 0.1   # width for smearing
	h = 0.05 # meshsize for the density of states
    x, y = dos(λ, σ, h)
    push!(X_W, x[1:100])
	push!(DoS_W, y[1:100])
end
plot(X_W, DoS_W, label = ["EcL = 1000, EcW = 5" "EcL = 1000, EcW = 10" "EcL = 1000, EcW = 50"], ylabel="DoS", lw = 3)



