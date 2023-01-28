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

# Generate EcL 
Wstep = 5
lengthW= 7
vecEcW = zeros(Int64,lengthW);
intialW = 5;
for i = 1 : lengthW
    vecEcW[i]=intialW + Wstep*(i-1)
end

g_fermi_beta001mu10(x) = x*(1+exp((x-10)*0.01))^(-1)
g_fermi_beta002mu10(x) = x*(1+exp((x-10)*0.02))^(-1)
g_fermi_beta01mu10(x) = x*(1+exp((x-10)*0.1))^(-1)

dosk_Fermi_beta001mu10_gamma001 = zeros(Float64,lengthW);
dosk_Fermi_beta002mu10_gamma001 = zeros(Float64,lengthW);
dosk_Fermi_beta01mu10_gamma001 = zeros(Float64,lengthW);

for i=1:lengthW
    EcL =5
    EcW =vecEcW[i]
	neigs=10
	nfx =4* floor( Int, (2.0 * EcL) / norm(B1[:,1]) )    # 4* Gmax11
	nfy = 4* floor( Int, (2.0 * EcL) / norm(B1[:,2]) )    # 4* Gmax12
    γ=0.01 
	# K-points sampling
    kpoint=[]
    h = 5
    nk = Int(2*EcW/h)
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
	""" β = 0.01, 0.02, 0.1"""
	gH_zero_beta001=zeros(Float64,Nk,size(λ,2));
	dosk_beta001=zeros(Float64,Nk)
	gH_zero_beta002=zeros(Float64,Nk,size(λ,2));
	dosk_beta002=zeros(Float64,Nk)
	gH_zero_beta01=zeros(Float64,Nk,size(λ,2));
	dosk_beta01=zeros(Float64,Nk)
	for k=1 : Nk
		for i=1:size(λ,2)
			gH_zero_beta001[k,i]=g_fermi_beta001mu10(λ[k,i])*norm(Ψ[k,zero_index,i])^2
			gH_zero_beta002[k,i]=g_fermi_beta002mu10(λ[k,i])*norm(Ψ[k,zero_index,i])^2
			gH_zero_beta01[k,i]=g_fermi_beta01mu10(λ[k,i])*norm(Ψ[k,zero_index,i])^2
		end
		dosk_beta001[k]=sum(gH_zero_beta001[k,:])
		dosk_beta002[k]=sum(gH_zero_beta002[k,:])
		dosk_beta01[k]=sum(gH_zero_beta01[k,:])
	end   
	dosk_Fermi_beta001mu10_gamma001[i]=sum(dosk_beta001)*h^2
	dosk_Fermi_beta002mu10_gamma001[i]=sum(dosk_beta002)*h^2
	dosk_Fermi_beta01mu10_gamma001[i]=sum(dosk_beta01)*h^2
end

# Exact solution

EcW = 60
# K-points sampling
kpoint=[]
h = 5
nk = Int(2*EcW/h)
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
""" β = 0.01, 0.02, 0.1"""
gH_zero_beta001_exact=zeros(Float64,Nk,size(λ,2));
dosk_beta001_exact=zeros(Float64,Nk)
gH_zero_beta002_exact=zeros(Float64,Nk,size(λ,2));
dosk_beta002_exact=zeros(Float64,Nk)
gH_zero_beta01_exact=zeros(Float64,Nk,size(λ,2));
dosk_beta01_exact=zeros(Float64,Nk)
for k=1 : Nk
	for i=1:size(λ,2)
		gH_zero_beta001_exact[k,i]=g_fermi_beta001mu10(λ[k,i])*norm(Ψ[k,zero_index,i])^2
		gH_zero_beta002_exact[k,i]=g_fermi_beta002mu10(λ[k,i])*norm(Ψ[k,zero_index,i])^2
		gH_zero_beta01_exact[k,i]=g_fermi_beta01mu10(λ[k,i])*norm(Ψ[k,zero_index,i])^2
	end
	dosk_beta001_exact[k]=sum(gH_zero_beta001_exact[k,:])
	dosk_beta002_exact[k]=sum(gH_zero_beta002_exact[k,:])
	dosk_beta01_exact[k]=sum(gH_zero_beta01_exact[k,:])
end   
dosk_Fermi_beta001mu10_gamma001_exact[i]=sum(dosk_beta001_exact)*h^2
dosk_Fermi_beta002mu10_gamma001_exact[i]=sum(dosk_beta002_exact)*h^2
dosk_Fermi_beta01mu10_gamma001_exact[i]=sum(dosk_beta01_exact)*h^2

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
