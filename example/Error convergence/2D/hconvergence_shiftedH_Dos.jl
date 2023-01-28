"""W = 30"""
using Plots
using LaTeXStrings

include("../../../src/struct.jl")
include("../../../src/2Dincommensurate_shiftedH.jl")

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

nkstep=5
lengthnk=10
vecnk=zeros(Int64,lengthnk);
intialnk=5;
for i=1:lengthnk 
    vecnk[i]=intialnk+nkstep*(i-1)
end

g_fermi_beta1mu10(x) = x*(1+exp((x-10)))^(-1)
g_fermi_beta01mu10(x) = x*(1+exp((x-10)*0.1))^(-1)

dosk_Fermi_beta1mu10_gamma005 = zeros(Float64,lengthnk);
dosk_Fermi_beta01mu10_gamma005 = zeros(Float64,lengthnk);

for i=1:lengthnk
    EcL = 5.
    EcW = 30.
    nk = vecnk[i]
    h = 2*EcW/nk 
    neigs = 10
    nfx = 4* floor( Int, (2.0 * EcL) / norm(B1[:,1]) )    # 4* Gmax11
    nfy= 4* floor( Int, (2.0 * EcL) / norm(B1[:,2]) )    # 4* Gmax12
    γ=0.05
    
    # K-points sampling
    kpoint=[]
    kpointx = range(-EcW, EcW-h, length=nk)
    for i=1:nk
        for j=1:nk
            j1= kpointx[i]
            j2= kpointx[j]
            push!(kpoint, [j1 j2]')
        end
    end
   Nk=size(kpoint,1)
   kpts = zeros(Float64,2,Nk)
   for i=1:Nk
       kpts[:,i] = kpoint[i][1:2]
   end  
   model = pwIncommensurate2D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftwx=nfx, n_fftwy=nfy, n_eigs=neigs,γ=γ)  

   # hamiltonian
   H, G, Gmn, R, Gmax11,Gmax12,Gmax21,Gmax22=hamiltonian2d(atoms, model);
   zero_index=Int64
   for i=1:size(G,1)
	   if G[i][:]==[0; 0; 0; 0]
		  zero_index=i
	   end
   end
   λ = zeros(Float64, Nk, size(H,2));
   Ψ = zeros(Float64, Nk, size(H,2), size(H,2));
   for k=1 : Nk
	   λ[k,:],Ψ[k,:,:] = eigen(H[k,:,:]); 
   end
   """ β = 1, 0.1"""
   gH_zero_beta1=zeros(Float64,Nk,size(λ,2));
   dosk_beta1=zeros(Float64,Nk)
   gH_zero_beta01=zeros(Float64,Nk,size(λ,2));
   dosk_beta01=zeros(Float64,Nk)
   for k=1 : Nk
	   for i=1:size(λ,2)
		   gH_zero_beta1[k,i]=g_fermi_beta1mu10(λ[k,i])*norm(Ψ[k,zero_index,i])^2
		   gH_zero_beta01[k,i]=g_fermi_beta01mu10(λ[k,i])*norm(Ψ[k,zero_index,i])^2
	   end
	   dosk_beta1[k]=sum(gH_zero_beta1[k,:])
	   dosk_beta01[k]=sum(gH_zero_beta01[k,:])
   end   
   dosk_Fermi_beta1mu10_gamma005[i]=sum(dosk_beta1)*h^2
   dosk_Fermi_beta01mu10_gamma005[i]=sum(dosk_beta01)*h^2
end

# Exact solution
nk= 100
h = 2*EcW/nk 
# K-points sampling
kpoint=[]
kpointx = range(-EcW, EcW-h, length=nk)
for i=1:nk
	for j=1:nk
		j1= kpointx[i]
		j2= kpointx[j]
		push!(kpoint, [j1 j2]')
	end
end
Nk=size(kpoint,1)
kpts = zeros(Float64,2,Nk)
for i=1:Nk
   kpts[:,i] = kpoint[i][1:2]
end  
model = pwIncommensurate2D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftwx=nfx, n_fftwy=nfy, n_eigs=neigs,γ=γ)  

# hamiltonian
H, G, Gmn, R, Gmax11,Gmax12,Gmax21,Gmax22=hamiltonian2d(atoms, model);
zero_index=Int64
for i=1:size(G,1)
   if G[i][:]==[0; 0; 0; 0]
	  zero_index=i
   end
end
λ = zeros(Float64, Nk, size(H,2));
Ψ = zeros(Float64, Nk, size(H,2), size(H,2));
for k=1 : Nk
   λ[k,:],Ψ[k,:,:] = eigen(H[k,:,:]); 
end
""" β = 1, 0.1"""
gH_zero_beta1_exact=zeros(Float64,Nk,size(λ,2));
dosk_beta1_exact=zeros(Float64,Nk)
gH_zero_beta01_exact=zeros(Float64,Nk,size(λ,2));
dosk_beta01_exact=zeros(Float64,Nk)
for k=1 : Nk
   for i=1:size(λ,2)
	   gH_zero_beta1_exact[k,i]=g_fermi_beta1mu10(λ[k,i])*norm(Ψ[k,zero_index,i])^2
	   gH_zero_beta01_exact[k,i]=g_fermi_beta01mu10(λ[k,i])*norm(Ψ[k,zero_index,i])^2
   end
   dosk_beta1_exact[k]=sum(gH_zero_beta1_exact[k,:])
   dosk_beta01_exact[k]=sum(gH_zero_beta01_exact[k,:])
end   
dosk_Fermi_beta1mu10_gamma005_exact=sum(dosk_beta1_exact)*h^2
dosk_Fermi_beta01mu10_gamma005_exact=sum(dosk_beta01_exact)*h^2

# Error for different β
herror_Fermi_beta1mu10_gamma005 = zeros(Float64,lengthnk);
herror_Fermi_beta01mu10_gamma005 = zeros(Float64,lengthnk);
for i = 1 : lengthnk
	herror_Fermi_beta1mu10_gamma005[i]=abs(dosk_Fermi_beta1mu10_gamma005_exact - dosk_Fermi_beta1mu10_gamma005[i])
	herror_Fermi_beta01mu10_gamma005[i]=abs(dosk_Fermi_beta01mu10_gamma005_exact  - dosk_Fermi_beta01mu10_gamma005[i])
end

# plot log - log 
vech = 2 * EcW * vecnk.^(-1)
vechinv = vech.^(-1)
herror = [herror_Fermi_beta1mu10_gamma005,herror_Fermi_beta01mu10_gamma005]
plot(vechinv, herror, yaxis=:log, label=[L"\beta=1" L"\beta=0.1"], xlabel=L"h^{-1}", ylabel= "Error", legend= :bottomleft, shape=[:circle :utriangle], markersize=6, lw=3, tickfontsize=15, legendfontsize=13, guidefontsize=15)