using Plots
using LaTeXStrings

include("../../src/struct.jl")
include("../../src/2Dincommensurate.jl")
include("../../src/dos.jl")

# set atom parameters 2d
X1 = [0.0;0.0]
X2 = [0.0;0.0]
s=1  
θ=π/10     # angle of rotation
L=2
A1=[1 0; 0 1]  
R1 =L*A1     #   sheet 1
A2=[cos(θ) cos(θ+π/2); sin(θ) sin(θ+π/2)] 
R2 =L*A2   # sheet 2 is rotated by θ=π/10
B1,B2=reciprocal(R1,R2)
atoms= TwoLayerIn2D(R1=R1, R2=R2,B1=B1,B2=B2,θ=θ, X1=X1, X2=X2)

# area of the reciprocal unit cell
RS1 = sqrt( norm(B1[:,1])^2 * norm(B1[:,2])^2 - dot(B1[:,1], B1[:,2])^2 ) 
RS2 = sqrt( norm(B2[:,1])^2 * norm(B2[:,2])^2 - dot(B2[:,1], B2[:,2])^2 )
Coe=RS1*RS2*π/4

# set the model parameters
EcL = 100.0
EcW = 5.0
neigs=10
nfx =4* floor( Int, (2.0 * EcL) / norm(B1[:,1]) )    # 4* Gmax11
nfy= 4* floor( Int, (2.0 * EcL) / norm(B1[:,2]) )    # 4* Gmax12
γ=0.01
kpts=[0.0 0.0;0.0 0.0]
model = pwIncommensurate2D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftwx=nfx, n_fftwy=nfy, n_eigs=neigs, γ=γ)  

H, G, Gmn, R, Gmax11,Gmax12,Gmax21,Gmax22=hamiltonian2d(atoms, model)
D, Ψ = eigen(Array(H));
λ = real.(D)
σ = 0.1   # width for smearing
h = 0.05 # meshsize for the density of states
x, y = dos(λ, σ, h)
# plot the first 100 points
plot(x[1:100], y[1:100], label = "EcL = 100, EcW = 5", ylabel="DoS", lw = 3)


# Dos convergence with respect to L
vecEcL = [50, 50, 80, 100]
EcW = 5.0

X_L = []
DoS_L = []
for i = 1 : length(vecEcL)
	EcL = vecEcL[i]
	model = pwIncommensurate2D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftwx=nfx, n_fftwy=nfy, n_eigs=neigs, γ=γ)    

	H, G, Gmn, R, Gmax11,Gmax12,Gmax21,Gmax22=hamiltonian2d(atoms, model)
	D, Ψ = eigen(Array(H));
	λ = real.(D)
	σ = 0.1   # width for smearing
	h = 0.05 # meshsize for the density of states
    x, y = dos(λ, σ, h)
    push!(X_L, x)
	push!(DoS_L, y)
end

plot(X_L, DoS_L, label = ["EcL = 50, EcW = 5" "EcL = 60, EcW = 5" "EcL = 80, EcW = 5" "EcL = 100, EcW = 5"], ylabel="DoS", lw = 3)


# Dos convergence with respect to W
EcL = 50 
vecEcW = [5, 8, 10];   # 30, 40,

X_W = []
DoS_W = []
for i = 1 : length(vecEcW)
	EcW = vecEcW[i]
	model = pwIncommensurate2D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftwx=nfx, n_fftwy=nfy, n_eigs=neigs, γ=γ) 

	H, G, Gmn, R, Gmax11,Gmax12,Gmax21,Gmax22=hamiltonian2d(atoms, model)
	D, Ψ = eigen(Array(H));
	λ = real.(D)
	σ = 0.1   # width for smearing
	h = 0.05 # meshsize for the density of states
    x, y = dos(λ, σ, h)
    push!(X_W, x[1:100])
	push!(DoS_W, y[1:100])
end
plot(X_W, DoS_W, label = ["EcL = 50, EcW = 5" "EcL = 50, EcW = 8" "EcL = 50, EcW = 10"], ylabel="DoS", lw = 3)

