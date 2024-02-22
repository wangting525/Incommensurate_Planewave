# Incommensurate_Planewave.jl 

`Incommensurate_Planewave.jl` is a Julia package designed for the computation of the density of states (DoS) and local density of states (LDoS) for 1D and 2D incommensurate systems. 

Including:

- Plane wave method.
- Hamiltonian Construction.
- Eigenpair Computation.
- DoS and LDoS Computation.

For an in-depth understanding of the algorithms employed, please refer to the article [Convergence of the planewave
approximations for quantum incommensurate systems. arXiv:2204.00994, 2022.].

Quick Strart 
***
Here is an example of how to perform basic calculations using `Incommensurate_Planewave.jl`:
```julia

using Plots
using LaTeXStrings

include("../../src/struct.jl")
include("../../src/2Dincommensurate.jl")
include("../../src/dos.jl")
```

2D incommensurate systems that are obtained by two periodic lattices together, in which one layer is rotated by an angle $\theta = \pi/10$ for the other.

```julia
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
```
Set the plane wave cutoff parameters.

```julia
# set the model parameters
EcL = 100.0
EcW = 5.0
neigs=10
nfx =4* floor( Int, (2.0 * EcL) / norm(B1[:,1]) )    # 4* Gmax11
nfy= 4* floor( Int, (2.0 * EcL) / norm(B1[:,2]) )    # 4* Gmax12
γ=0.01
kpts=[0.0 0.0;0.0 0.0]
model = pwIncommensurate2D(EcL=EcL, EcW=EcW, kpts=kpts, n_fftwx=nfx, n_fftwy=nfy, n_eigs=neigs, γ=γ)  
```

Generate Hamiltonian matrix and solve eigenpairs. 
```julia
H, G, Gmn, R, Gmax11,Gmax12,Gmax21,Gmax22=hamiltonian2d(atoms, model)
D, Ψ = eigen(Array(H));
```

Compute the plot the density of states.
```julia
λ = real.(D)
σ = 0.1   # width for smearing
h = 0.05 # meshsize for the density of states
x, y = dos(λ, σ, h)
# plot the first 100 points
plot(x[1:100], y[1:100]./EcL^2, label = "EcL = 100, EcW = 5", ylabel="DoS", lw = 3)
```
