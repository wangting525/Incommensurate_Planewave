"""W = 5"""

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

Nkstep=5
lengthNk=20
vecNk=zeros(Int64,lengthNk);
intialNk=5;
for i=1:lengthNk  
    vecNk[i]=intialNk+Nkstep*(i-1)
end

g_fermi_beta1mu1(x) = x*(1+exp((x-1)*1))^(-1)
g_fermi_beta01mu1(x) = x*(1+exp((x-1)*0.1))^(-1)

dosk_Fermi_beta1mu1_gamma001 = zeros(Float64,lengthNk);
dosk_Fermi_beta01mu1_gamma001 = zeros(Float64,lengthNk);

for i=1:lengthNk
    EcL =25
    EcW =5
    neigs=10
    nf= 2 * floor( Int, 2*EcL *L1/π )  
    γ=0.01 
	nk = vecNk[i]
	h = 2 * EcW / nk
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
	""" β = 1, 0.1"""
	gH_zero_beta1=zeros(Float64,nk,size(λ,2));
	dosk_beta1=zeros(Float64,nk)
	gH_zero_beta01=zeros(Float64,nk,size(λ,2));
	dosk_beta01=zeros(Float64,nk)
	for k=1 : nk
		for i=1:size(λ,2)
			gH_zero_beta1[k,i]=g_fermi_beta1mu1(λ[k,i])*norm(Ψ[k,zero_index,i])^2
			gH_zero_beta01[k,i]=g_fermi_beta01mu1(λ[k,i])*norm(Ψ[k,zero_index,i])^2
		end
		dosk_beta1[k]=sum(gH_zero_beta1[k,:])
		dosk_beta01[k]=sum(gH_zero_beta01[k,:])
	end 
	dosk_Fermi_beta1mu1_gamma001[i]=sum(dosk_beta1)*h
	dosk_Fermi_beta01mu1_gamma001[i]=sum(dosk_beta01)*h
end

# Exact solution
nk = 200
h = 2 * EcW / nk
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
""" β = 1, 0.1"""
gH_zero_beta1_exact=zeros(Float64,nk,size(λ,2));
dosk_beta1_exact=zeros(Float64,nk)
gH_zero_beta01_exact=zeros(Float64,nk,size(λ,2));
dosk_beta01_exact=zeros(Float64,nk)
for k=1 : nk
	for i=1:size(λ,2)
		gH_zero_beta1_exact[k,i]=g_fermi_beta1mu1(λ[k,i])*norm(Ψ[k,zero_index,i])^2
		gH_zero_beta01_exact[k,i]=g_fermi_beta01mu1(λ[k,i])*norm(Ψ[k,zero_index,i])^2
	end
	dosk_beta1_exact[k]=sum(gH_zero_beta1_exact[k,:])
	dosk_beta01_exact[k]=sum(gH_zero_beta01_exact[k,:])
end 

dosk_Fermi_beta1mu1_gamma001_exact=sum(dosk_beta1_exact)*h
dosk_Fermi_beta01mu1_gamma001_exact=sum(dosk_beta01_exact)*h

# Error for different β
herror_Fermi_beta1mu1_gamma001 = zeros(Float64,lengthNk);
herror_Fermi_beta01mu1_gamma001 = zeros(Float64,lengthNk);
for i = 1 : lengthNk
	herror_Fermi_beta1mu1_gamma001[i]=abs(dosk_Fermi_beta1mu1_gamma001_exact - dosk_Fermi_beta1mu1_gamma001[i])
	herror_Fermi_beta01mu1_gamma001[i]=abs(dosk_Fermi_beta01mu1_gamma001_exact  - dosk_Fermi_beta01mu1_gamma001[i])
end

# plot log - log 
vech = 2 * EcW * vecNk.^(-1)
vechinv = vech.^(-1)
herror = [herror_Fermi_beta1mu1_gamma001,herror_Fermi_beta01mu1_gamma001]
plot(vechinv, herror, xaxis=:log, yaxis=:log, label=[L"\beta=1" L"\beta=0.1"], xlabel=L"h^{-1}", ylabel= "Error", legend= :bottomleft, shape=[:circle :utriangle], markersize=6, lw=3, tickfontsize=15, legendfontsize=15, guidefontsize=15)
plot!(vechinv, 20*vech, xaxis=:log, yaxis=:log, label= L"\mathcal{O}(L^{-1})", ls=:dash, lw=3, c=:black)