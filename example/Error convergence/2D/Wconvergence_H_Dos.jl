using Plots
using LaTeXStrings

include("../../../src/struct.jl")
include("../../../src/2Dincommensurate.jl")

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
	EcL = 20
	EcW = vecEcW[i]
	neigs=10
	nfx =4* floor( Int, (2.0 * EcL) / norm(B1[:,1]) )    # 4* Gmax11
	nfy= 4* floor( Int, (2.0 * EcL) / norm(B1[:,2]) )    # 4* Gmax12
	γ=0.001
	kpts=[0.0 0.0;0.0 0.0]
	model = pwIncommensurate2D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftwx=nfx, n_fftwy=nfy, n_eigs=neigs, γ=γ) 
    # hamiltonian
    H, G, Gmn, R, Gmax11,Gmax12,Gmax21,Gmax22=hamiltonian2d(atoms, model)
    λ = zeros(Float64, size(H,1));
    D, Ψ = eigen(Array(H));
    d_fermi_beta1mu10_EcW = zeros(Float64,length(D))
    for n = 1 : length(D)
        d_fermi_beta1mu10_EcW[n] = g(D[n]) 
    end
    dos_fermi_beta1mu10_gamma0001[i] = Coe*(1/EcL^2)*sum(d_fermi_beta1mu10_EcW)
end

# Exact solution
EcW = 50
model = pwIncommensurate2D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftwx=nfx, n_fftwy=nfy, n_eigs=neigs, γ=γ) 
# hamiltonian
H, G, Gmn, R, Gmax11,Gmax12,Gmax21,Gmax22=hamiltonian2d(atoms, model)
λ = zeros(Float64, size(H,1));
D, Ψ = eigen(Array(H));
d_fermi_beta1mu10_EcW_exact = zeros(Float64,length(D))
for n = 1 : length(D)
	d_fermi_beta1mu10_EcW_exact[n] = g(D[n]) 
end
dos_fermi_beta1mu10_gamma0001_exact = Coe*(1/EcL^2)*sum(d_fermi_beta1mu10_EcW_exact)

"""γ = 0.01"""
dos_fermi_beta1mu10_gamma001 = zeros(Float64,lengthW);
for i=1:lengthW
	EcL = 20
	EcW = vecEcW[i]
	neigs=10
	nfx =4* floor( Int, (2.0 * EcL) / norm(B1[:,1]) )    # 4* Gmax11
	nfy= 4* floor( Int, (2.0 * EcL) / norm(B1[:,2]) )    # 4* Gmax12
	γ=0.01
	kpts=[0.0 0.0;0.0 0.0]
	model = pwIncommensurate2D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftwx=nfx, n_fftwy=nfy, n_eigs=neigs, γ=γ) 
    # hamiltonian
    H, G, Gmn, R, Gmax11,Gmax12,Gmax21,Gmax22=hamiltonian2d(atoms, model)
    λ = zeros(Float64, size(H,1));
    D, Ψ = eigen(Array(H));
    d_fermi_beta1mu10_EcW = zeros(Float64,length(D))
    for n = 1 : length(D)
        d_fermi_beta1mu10_EcW[n] = g(D[n]) 
    end
    dos_fermi_beta1mu10_gamma001[i] = Coe*(1/EcL^2)*sum(d_fermi_beta1mu10_EcW)
end

# Exact solution
EcW = 50
model = pwIncommensurate2D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftwx=nfx, n_fftwy=nfy, n_eigs=neigs, γ=γ) 
# hamiltonian
H, G, Gmn, R, Gmax11,Gmax12,Gmax21,Gmax22=hamiltonian2d(atoms, model)
λ = zeros(Float64, size(H,1));
D, Ψ = eigen(Array(H));
d_fermi_beta1mu10_EcW_exact = zeros(Float64,length(D))
for n = 1 : length(D)
	d_fermi_beta1mu10_EcW_exact[n] = g(D[n]) 
end
dos_fermi_beta1mu10_gamma001_exact = Coe*(1/EcL^2)*sum(d_fermi_beta1mu10_EcW_exact)

"""γ = 0.1"""
dos_fermi_beta1mu10_gamma01 = zeros(Float64,lengthW);
for i=1:lengthW
	EcL = 20
	EcW = vecEcW[i]
	neigs=10
	nfx =4* floor( Int, (2.0 * EcL) / norm(B1[:,1]) )    # 4* Gmax11
	nfy= 4* floor( Int, (2.0 * EcL) / norm(B1[:,2]) )    # 4* Gmax12
	γ=0.1
	kpts=[0.0 0.0;0.0 0.0]
	model = pwIncommensurate2D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftwx=nfx, n_fftwy=nfy, n_eigs=neigs, γ=γ) 
    # hamiltonian
    H, G, Gmn, R, Gmax11,Gmax12,Gmax21,Gmax22=hamiltonian2d(atoms, model)
    λ = zeros(Float64, size(H,1));
    D, Ψ = eigen(Array(H));
    d_fermi_beta1mu10_EcW = zeros(Float64,length(D))
    for n = 1 : length(D)
        d_fermi_beta1mu10_EcW[n] = g(D[n]) 
    end
    dos_fermi_beta1mu10_gamma01[i] = Coe*(1/EcL^2)*sum(d_fermi_beta1mu10_EcW)
end

# Exact solution
EcW = 50
model = pwIncommensurate2D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftwx=nfx, n_fftwy=nfy, n_eigs=neigs, γ=γ) 
# hamiltonian
H, G, Gmn, R, Gmax11,Gmax12,Gmax21,Gmax22=hamiltonian2d(atoms, model)
λ = zeros(Float64, size(H,1));
D, Ψ = eigen(Array(H));
d_fermi_beta1mu10_EcW_exact = zeros(Float64,length(D))
for n = 1 : length(D)
	d_fermi_beta1mu10_EcW_exact[n] = g(D[n]) 
end
dos_fermi_beta1mu10_gamma01_exact = Coe*(1/EcL^2)*sum(d_fermi_beta1mu10_EcW_exact)

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