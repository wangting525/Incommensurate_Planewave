using Plots
using LaTeXStrings

include("../../../src/struct.jl")
include("../../../src/1Dincommensurate_shiftedH.jl")

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

g_fermi_beta001mu55(x) = x*(1+exp((x-55)*0.01))^(-1)
g_fermi_beta002mu55(x) = x*(1+exp((x-55)*0.02))^(-1)
g_fermi_beta01mu55(x) = x*(1+exp((x-55)*0.1))^(-1)

# Generate EcW 
Wstep = 5
lengthW = 7
vecEcW = zeros(Int64,lengthW);
intialW = 5;
for i=1:lengthW
    vecEcW[i]=intialW+Wstep*(i-1)
end

dosk_Fermi_beta001mu10_gamma001 = zeros(Float64,lengthW);
dosk_Fermi_beta002mu10_gamma001 = zeros(Float64,lengthW);
dosk_Fermi_beta01mu10_gamma001 = zeros(Float64,lengthW);
for i=1:lengthW
    EcL =20
    EcW =vecEcW[i]
    neigs=10
    nf= 2 * floor( Int, 2*EcL *L1/π )  
    γ=0.01 
	h = 0.5
	nk = floor(Int, 2 * EcW / h)
    kpoint = range(-EcW, EcW-h, length=nk)
    kpts = [kpoint[i] for i=1:length(kpoint)]
    model = pwIncommensurate1D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ)
    # hamiltonian
    H, G, R=hamiltonian(atoms, model);
	zero_index=Int64
	for i=1:length(G)
		if G[i][1]==0
		   if G[i][2]==0
			  zero_index=i
		   end
		end
	end
	λ = zeros(Float64, nk, size(H,2));
	Ψ = zeros(Float64, nk, size(H,2), size(H,2));
	for k=1 : nk
		λ[k,:],Ψ[k,:,:] = eigen(H[k,:,:]); 
	end
	""" β = 0.01, 0.02, 0.1"""
	gH_zero_beta001=zeros(Float64,nk,size(λ,2));
	dosk_beta001=zeros(Float64,nk)
	gH_zero_beta002=zeros(Float64,nk,size(λ,2));
	dosk_beta002=zeros(Float64,nk)
	gH_zero_beta01=zeros(Float64,nk,size(λ,2));
	dosk_beta01=zeros(Float64,nk)
	for k=1 : nk
		for i=1:size(λ,2)
			gH_zero_beta001[k,i]=g_fermi_beta001mu55(λ[k,i])*norm(Ψ[k,zero_index,i])^2
			gH_zero_beta002[k,i]=g_fermi_beta002mu55(λ[k,i])*norm(Ψ[k,zero_index,i])^2
			gH_zero_beta01[k,i]=g_fermi_beta01mu55(λ[k,i])*norm(Ψ[k,zero_index,i])^2
		end
		dosk_beta001[k]=sum(gH_zero_beta001[k,:])
		dosk_beta002[k]=sum(gH_zero_beta002[k,:])
		dosk_beta01[k]=sum(gH_zero_beta01[k,:])
	end   
	dosk_Fermi_beta001mu10_gamma001[i]=sum(dosk_beta001)*h
	dosk_Fermi_beta002mu10_gamma001[i]=sum(dosk_beta002)*h
	dosk_Fermi_beta01mu10_gamma001[i]=sum(dosk_beta01)*h
end

# Exact solution
EcW =100
γ=0.01 
h = 0.5
nk = floor(Int, 2 * EcW / h)
kpoint = range(-EcW, EcW-h, length=nk)
kpts = [kpoint[i] for i=1:length(kpoint)]
model = pwIncommensurate1D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ)
# hamiltonian
H, G, R=hamiltonian(atoms, model);
zero_index=Int64
for i=1:length(G)
	if G[i][1]==0
	   if G[i][2]==0
		  zero_index=i
	   end
	end
end
λ = zeros(Float64, nk, size(H,2));
Ψ = zeros(Float64, nk, size(H,2), size(H,2));
for k=1 : nk
	λ[k,:],Ψ[k,:,:] = eigen(H[k,:,:]); 
end
""" β = 0.01, 0.02, 0.1"""
gH_zero_beta001_exact=zeros(Float64,nk,size(λ,2));
dosk_beta001_exact=zeros(Float64,nk)
gH_zero_beta002_exact=zeros(Float64,nk,size(λ,2));
dosk_beta002_exact=zeros(Float64,nk)
gH_zero_beta01_exact=zeros(Float64,nk,size(λ,2));
dosk_beta01_exact=zeros(Float64,nk)
for k=1 : nk
	for i=1:size(λ,2)
		gH_zero_beta001_exact[k,i]=g_fermi_beta001mu55(λ[k,i])*norm(Ψ[k,zero_index,i])^2
		gH_zero_beta002_exact[k,i]=g_fermi_beta002mu55(λ[k,i])*norm(Ψ[k,zero_index,i])^2
		gH_zero_beta01_exact[k,i]=g_fermi_beta01mu55(λ[k,i])*norm(Ψ[k,zero_index,i])^2
	end
	dosk_beta001_exact[k]=sum(gH_zero_beta001[k,:])
	dosk_beta002_exact[k]=sum(gH_zero_beta002[k,:])
	dosk_beta01_exact[k]=sum(gH_zero_beta01[k,:])
end   
dosk_Fermi_beta001mu10_gamma001_exact=sum(dosk_beta001)*h
dosk_Fermi_beta002mu10_gamma001_exact=sum(dosk_beta002)*h
dosk_Fermi_beta01mu10_gamma001_exact=sum(dosk_beta01)*h

# Error for different β
Werror_Fermi_beta001mu10_gamma001 = zeros(Float64,lengthW);
Werror_Fermi_beta002mu10_gamma001 = zeros(Float64,lengthW);
Werror_Fermi_beta01mu10_gamma001 = zeros(Float64,lengthW);
for i = 1 : lengthW
	Werror_Fermi_beta001mu10_gamma001[i]=abs(dosk_Fermi_beta001mu10_gamma001_exact - dosk_Fermi_beta001mu10_gamma001[i])
	Werror_Fermi_beta002mu10_gamma001[i]=abs(dosk_Fermi_beta002mu10_gamma001_exact  - dosk_Fermi_beta002mu10_gamma001[i])
	Werror_Fermi_beta001mu10_gamma001[i]=abs(dosk_Fermi_beta01mu10_gamma001_exact  - dosk_Fermi_beta01mu10_gamma001[i])
end

Werror = [Werror_Fermi_beta001mu10_gamma001,Werror_Fermi_beta002mu10_gamma001,Werror_Fermi_beta01mu10_gamma001]
plot(vecEcW, Werror, yaxis=:log, label=[L"\beta=0.001" L"\beta=0.02" L"\beta=0.1"],  
     xlabel="W", ylabel= "Error", legend= :bottomleft,
 shape=[:circle :utriangle :rect], markersize=6, lw=3, tickfontsize=15, legendfontsize=15, guidefontsize=15)
