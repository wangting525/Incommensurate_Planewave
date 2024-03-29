# Incommensurate_Planewave.jl 

`Incommensurate_Planewave.jl` is a Julia package designed for the computation of the density of states (DoS) and local density of states (LDoS) for 1D and 2D incommensurate systems. 

Including:

- Plane Wave Approximations.
- Hamiltonian Construction.
- Eigenpair Computation.
- DoS and LDoS Computation.

For detailed numerical schemes, refer to the article below:

T. Wang, H. Chen, A. Zhou, Y. Zhou and D. Massatt, <a href="https://arxiv.org/abs/2204.00994" style="color: blue;">Convergence of the planewave
approximations for quantum incommensurate systems</a>, arXiv:2204.00994, 2022.

## Quick Start 

Example of basic calculations  using `Incommensurate_Planewave.jl`:
```julia

using Plots
using LaTeXStrings

include("../../src/struct.jl")
include("../../src/2Dincommensurate.jl")
include("../../src/dos.jl")
```

2D incommensurate systems that are obtained by two periodic lattices together, in which one layer is rotated by an angle $\theta$ for the other.

<img src="https://github.com/wangting525/Incommensurate_Planewave/blob/master/figures/atomstructure.png" width="300" alt="Incommensurate system">


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
H = hamiltonian2d(atoms, model)
D, Ψ = eigen(Array(H));
```

Compute and plot the density of states.
```julia
λ = real.(D)
σ = 0.1   # width for smearing
h = 0.05 # meshsize for the density of states
x, y = dos(λ, σ, h)
# plot the first 100 points
plot(x[1:100], y[1:100]./EcL^2, label = "EcL = 100, EcW = 5", ylabel="DoS", lw = 3)
```
<img src="https://github.com/wangting525/Incommensurate_Planewave/blob/master/figures/1Ddos.png" width="500" alt="Incommensurate system">
