using StaticArrays
using SparseArrays
using LinearAlgebra
using Arpack
using FFTW

"""`Hamiltonian`"""
function hamiltonian2d(atoms::TwoLayerIn2D, model::pwIncommensurate2D)
    R1 = atoms.R1
    R2 = atoms.R2
    B1 = atoms.B1
    B2 = atoms.B2

    kpts = model.kpts
    EcL = model.EcL
    EcW = model.EcW
    γ=model.γ

    Gmax11 = floor( Int, max(EcL,EcW) / norm(B1[:,1]) )
    Gmax12 = floor( Int, max(EcL,EcW) / norm(B1[:,2]) )
    Gmax21 = floor( Int, max(EcL,EcW) / norm(B2[:,1]) )
    Gmax22 = floor( Int, max(EcL,EcW) / norm(B2[:,2]) )

    # area of the primative unit cell
    S1 = sqrt( norm(R1[:,1])^2 * norm(R1[:,2])^2 - dot(R1[:,1], R1[:,2])^2 )
    S2 = sqrt( norm(R2[:,1])^2 * norm(R2[:,2])^2 - dot(R2[:,1], R2[:,2])^2 )

    #area of the reciprocal unit cell
    RS1 = sqrt( norm(B1[:,1])^2 * norm(B1[:,2])^2 - dot(B1[:,1], B1[:,2])^2 )
    RS2 = sqrt( norm(B2[:,1])^2 * norm(B2[:,2])^2 - dot(B2[:,1], B2[:,2])^2 )

    # count the degrees of freedom and label the basis functions
    # compute npw
    G=[]
    Gmn=[]
    npw = 0
    for j11 = -Gmax11 : Gmax11, j12 = -Gmax12 : Gmax12
        for j21 = -Gmax21 : Gmax21, j22 = -Gmax22 : Gmax22
            G1m = j11 * B1[:,1] + j12 * B1[:,2];    #  G1m=[G1_m1;G1_m2]
            G2n = j21 * B2[:,1] + j22 * B2[:,2];
            if norm( G1m + G2n ) < EcW
                if norm( G1m - G2n) < EcL
                   npw = npw + 1
                   push!(G, [j11 j12 j21 j22]')
                   push!(Gmn,[G1m[1,1] G1m[2,1] G2n[1,1] G2n[2,1]]')
               end
           end
        #end
        end
    end
    println(" EcutL = ", EcL, "; EcutW = ", EcW, "; DOF = ", npw)
    # assemble the matrix
    nk = size(kpts,2)
    H = zeros(Complex{Float64}, nk, npw, npw)
    R = zeros(Float64,2,npw)
    #diffG=zeros(Int64,1,npw*npw)
 for k = 1 : nk
     println("assembling H[k] with kpt = ", kpts[:,k])

    for G1 = 1 : npw, G2 = 1 : npw
        diffG=G[G1][:]-G[G2][:]
        # sample k points
        #NormR[1,G1]=norm(Gmn[:,G1] + kpt[:])^2 +
        #2*dot( Gmn[1:2,G1] + kpt[:,1] , Gmn[3:4,G1] + kpt[:,2] )
        H[k,G1, G2] = ( G1==G2 ?
             #0.5 * (norm(Gmn[1:2,G1] + Gmn[3:4,G1] + kpts[:,k])^2 )
             #+1 +1
             0.5 *( norm(kpts[:,k])^2 + norm(Gmn[G1][:])^2) + dot(kpts[:,k] + Gmn[G1][1:2], Gmn[G1][3:4]) + 1 + 1
            # + 2*dot( Gmn[1:2,G1] + kpt[:,:], Gmn[3:4,G1] + kpt[:,2]))+1+1 :
            : ( G[G1][3:4] == G[G2][3:4] ?
                     exp( -γ * norm(B1 * diffG[1:2])^2 ) : 0  )
                     #vfft1[ mod((G[1,G1]-G[1,G2]),nfx)+1, mod((G[2,G1]-G[2,G2]),nfy)+1 ] :  0 )
                        # * exp( im * dot(G[1:2,G2]-G[1:2,G1], Ratom) )
            + ( G[G1][1:2] == G[G2][1:2] ?
                     exp( -γ * norm(B2 * diffG[3:4])^2 ) : 0  )
                    #vfft2[ mod((G[3,G1]-G[3,G2]),nfx)+1, mod((G[4,G1]-G[4,G2]),nfy)+1 ] : 0 )
                        # * exp( im * dot(G[3:4,G2]-G[3:4,G1], Ratom) )
             )
            R[:,G1]=Gmn[G1][1:2]+Gmn[G1][3:4]     # G1m+G2n
     end
  end
    return H, G, Gmn, R ,Gmax11, Gmax12, Gmax21, Gmax22, S1, S2, RS1, RS2  #,vfft1,vfft2
end
