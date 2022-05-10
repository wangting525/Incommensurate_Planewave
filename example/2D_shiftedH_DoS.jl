using Plots
using LaTeXStrings

include("../src/struct.jl")
include("../src/2Dincommensurate_shiftedH.jl")

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
v(x,y)=s*((cos((π/L)*x))^2)*((cos((π/L)*y))^2)
v1(x,y)=v((A1*[x;y])[1],(A1*[x;y])[2])
v2(x,y)=v((A2*[x;y])[1],(A2*[x;y])[2])
#v1(x,y)=s1*(cos(π*(A1*[x;y])[1]))^2*(cos(π*(A1*[x;y])[2]))^2  
#v2(x,y)=s2*(cos(π*(A2*[x;y])[1]))^2*(cos(π*(A2*[x;y])[2]))^2 
B1,B2=reciprocal(R1,R2)
atoms= TwoLayerIn2D(R1=R1, R2=R2,B1=B1,B2=B2,v1=v1, v2=v2,θ=θ, X1=X1, X2=X2)

#set the model parameters
EcL = 10.0
EcW = 5.0
neigs=10
nfx =4* floor( Int, (2.0 * EcL) / norm(B1[:,1]) )    # 4* Gmax11
nfy= 4* floor( Int, (2.0 * EcL) / norm(B1[:,2]) )    # 4* Gmax12
γ=0.05
#K-points uniform sampling
nk=5
h=2*EcW/nk
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

# test function
g_fermi_beta01mu10(x) = x*(1+exp((x-10)*0.1))^(-1)
g_fermi_beta1mu10(x) = x*(1+exp((x-10)*1))^(-1)
g_fermi_beta10mu10(x) = x*(1+exp((x-10)*10))^(-1)
g(x)=g_fermi_beta1mu10(x)

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

# DoS for test function with β=1,μ=10 
gH_zero=zeros(Float64,Nk,size(λ,2));
dos_nk=zeros(Float64,Nk)
for k=1 : Nk
	for i=1:size(λ,2)
		gH_zero[k,i]=g(λ[k,i])*norm(Ψ[k,zero_index,i])^2
	end
	dos_nk[k]=sum(gH_zero[k,:])
end   
dos=sum(dos_nk)*h^2 