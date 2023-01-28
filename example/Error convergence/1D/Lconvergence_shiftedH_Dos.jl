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

# Generate EcL 
Lstep = 5
lengthL= 7
vecEcL = zeros(Int64,lengthL);
intialL = 5;
for i = 1 : lengthL
    vecEcL[i]=intialL + Lstep*(i-1)
end

g_fermi_beta10mu2(x) = x*(1+exp((x-2)*10))^(-1)
g_fermi_beta1mu2(x) = x*(1+exp((x-2)*1))^(-1)
g_fermi_beta01mu2(x) = x*(1+exp((x-2)*0.1))^(-1)

dosk_Fermi_beta10mu2_gamma001 = zeros(Float64,lengthL);
dosk_Fermi_beta1mu2_gamma001 = zeros(Float64,lengthL);
dosk_Fermi_beta01mu2_gamma001 = zeros(Float64,lengthL);
for i = 1 : lengthL
    EcL = vecEcL[i]
	EcW = 5
    neigs = 10
    nf = 2 * floor( Int, 2*EcL *L1/π ) 
    γ = 0.01 
	nk = 20
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
	""" β = 10, 1, 0.1"""
	gH_zero_beta10=zeros(Float64,nk,size(λ,2));
	dosk_beta10=zeros(Float64,nk)
	gH_zero_beta1=zeros(Float64,nk,size(λ,2));
	dosk_beta1=zeros(Float64,nk)
	gH_zero_beta01=zeros(Float64,nk,size(λ,2));
	dosk_beta01=zeros(Float64,nk)
	for k=1 : nk
		for i=1:size(λ,2)
			gH_zero_beta10[k,i]=g_fermi_beta10mu2(λ[k,i])*norm(Ψ[k,zero_index,i])^2
			gH_zero_beta1[k,i]=g_fermi_beta1mu2(λ[k,i])*norm(Ψ[k,zero_index,i])^2
			gH_zero_beta01[k,i]=g_fermi_beta01mu2(λ[k,i])*norm(Ψ[k,zero_index,i])^2
		end
		dosk_beta10[k]=sum(gH_zero_beta10[k,:])
		dosk_beta1[k]=sum(gH_zero_beta1[k,:])
		dosk_beta01[k]=sum(gH_zero_beta01[k,:])
	end   
	dosk_Fermi_beta10mu2_gamma001[i]=sum(dosk_beta10)*h
	dosk_Fermi_beta1mu2_gamma001[i]=sum(dosk_beta1)*h
	dosk_Fermi_beta01mu2_gamma001[i]=sum(dosk_beta01)*h
end

# Exact solution
EcL = 100
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
""" β = 10, 1, 0.1"""
gH_zero_beta10_exact=zeros(Float64,nk,size(λ,2));
dosk_beta10_exact=zeros(Float64,nk)
gH_zero_beta1_exact=zeros(Float64,nk,size(λ,2));
dosk_beta1_exact=zeros(Float64,nk)
gH_zero_beta01_exact=zeros(Float64,nk,size(λ,2));
dosk_beta01_exact=zeros(Float64,nk)
for k=1 : nk
	for i=1:size(λ,2)
		gH_zero_beta10_exact[k,i]=g_fermi_beta10mu2(λ[k,i])*norm(Ψ[k,zero_index,i])^2
		gH_zero_beta1_exact[k,i]=g_fermi_beta1mu2(λ[k,i])*norm(Ψ[k,zero_index,i])^2
		gH_zero_beta01_exact[k,i]=g_fermi_beta01mu2(λ[k,i])*norm(Ψ[k,zero_index,i])^2
	end
	dosk_beta10_exact[k]=sum(gH_zero_beta10[k,:])
	dosk_beta1_exact[k]=sum(gH_zero_beta1[k,:])
	dosk_beta01_exact[k]=sum(gH_zero_beta01[k,:])
end   
dosk_Fermi_beta10mu2_gamma001_exact=sum(dosk_beta10)*h
dosk_Fermi_beta1mu2_gamma001_exact=sum(dosk_beta1)*h
dosk_Fermi_beta01mu2_gamma001_exact=sum(dosk_beta01)*h

# Error for different β
Lerror_Fermi_beta10mu2_gamma001 = zeros(Float64,lengthL);
Lerror_Fermi_beta1mu2_gamma001 = zeros(Float64,lengthL);
Lerror_Fermi_beta01mu2_gamma001 = zeros(Float64,lengthL);
for i = 1 : lengthL
	Lerror_Fermi_beta10mu2_gamma001[i]=abs(dosk_Fermi_beta10mu2_gamma001_exact - dosk_Fermi_beta10mu2_gamma001[i])
	Lerror_Fermi_beta1mu2_gamma001[i]=abs(dosk_Fermi_beta1mu2_gamma001_exact  - dosk_Fermi_beta1mu2_gamma001[i])
	Lerror_Fermi_beta01mu2_gamma001[i]=abs(dosk_Fermi_beta01mu2_gamma001_exact  - dosk_Fermi_beta01mu2_gamma001[i])
end

Lerror = [Lerror_Fermi_beta10mu2_gamma001,Lerror_Fermi_beta1mu2_gamma001,Lerror_Fermi_beta01mu2_gamma001]
plot(vecEcL, Lerror, yaxis=:log, label=[L"\beta=10" L"\beta=1" L"\beta=0.1"],  
     xlabel="L", ylabel= "Error", legend= :bottomleft,
 shape=[:circle :utriangle :rect], markersize=6, lw=3, tickfontsize=15, legendfontsize=15, guidefontsize=15)
