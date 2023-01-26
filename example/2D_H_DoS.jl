using Plots
using LaTeXStrings

include("../src/struct.jl")
include("../src/2Dincommensurate.jl")

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
EcL = 20.0
EcW = 5.0
neigs=10
nfx =4* floor( Int, (2.0 * EcL) / norm(B1[:,1]) )    # 4* Gmax11
nfy= 4* floor( Int, (2.0 * EcL) / norm(B1[:,2]) )    # 4* Gmax12
γ=0.01
kpts=[0.0 0.0;0.0 0.0]
model = pwIncommensurate2D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftwx=nfx, n_fftwy=nfy, n_eigs=neigs, γ=γ)    

# test function
g_fermi_beta01mu10(x) = x*(1+exp((x-10)*0.1))^(-1)
g_fermi_beta1mu10(x) = x*(1+exp((x-10)*1))^(-1)
g_fermi_beta10mu10(x) = x*(1+exp((x-10)*10))^(-1)
g(x)=g_fermi_beta1mu10(x)

# hamiltonian
H, G, Gmn, R, Gmax11,Gmax12,Gmax21,Gmax22=hamiltonian2d(atoms, model)

λ = zeros(Float64, size(H,1));
D, Ψ = eigen(Array(H));

 # DoS for test function with β=1,μ=10 
 d=zeros(Float64,length(D))
 for i=1:length(D)
	 d[i]=g(D[i])
 end 
 dos=Coe*(1/EcL^2)*sum(d)