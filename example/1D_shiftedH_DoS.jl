using Plots
using LaTeXStrings

include("../src/struct.jl")
include("../src/1Dincommensurate_shiftedH.jl")

# set atom parameters
X1 = [0.0]
X2 = [0.0]
z = 0.0
s1= 1 
s2= 1 
L1= √5-1 
L2= 2
v1(x)=0.5*x^2
v2(x)=0.5*x^2 
atoms = BiLayer1D(L1=L1, L2=L2, X1=X1, X2=X2, z=z, v1=v1, v2=v2)

# set the model parameters 
EcL =20.0
EcW =10.0
nk = 20
neigs=10
nf= 2 * floor( Int, 2*EcL *L1/π )  
γ=0.01  
h = 2*EcW/(nk-1)
kpoint = range(-EcW, EcW, length=nk)
kpts = [kpoint[i] for i=1:length(kpoint)]
model = pwIncommensurate1D_LW(EcL=EcL, EcW=EcW, kpts=kpts, n_fftw=nf, n_eigs=neigs, γ=γ);

# test function
g_fermi_beta01mu10(x) = x*(1+exp((x-10)*0.1))^(-1)
g_fermi_beta1mu10(x) = x*(1+exp((x-10)*1))^(-1)
g_fermi_beta10mu10(x) = x*(1+exp((x-10)*10))^(-1)

g(x)=g_fermi_beta1mu10(x)

# hamiltonian
H ,G ,R = hamiltonian(atoms, model);

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
 
#Dos for test function with β=1,μ=10 
gH_zero=zeros(Float64,nk,size(λ,2));
dos_nk=zeros(Float64,nk)
for k=1 : nk
    for i=1:size(λ,2)
        gH_zero[k,i]=g(λ[k,i])*norm(Ψ[k,zero_index,i])^2
    end
    dos_nk[k]=sum(gH_zero[k,:])
end   
dos=sum(dos_nk)*h
