! Geometry block from utils/Main_BandStructure.f90, FORGE dev
! acbddd097c689142a38e5b0e99ac5fd674c0904a; in-plane assignment corrected.
aMoire = 3.0_dp*ntheta**2 + 3.0_dp*ntheta + 1.0_dp
g1  =  (4.0_dp*pi/3.0_dp)/aMoire*(real(3*ntheta+1,dp)*a1+a2)
g12 =  (4.0_dp*pi/3.0_dp)/aMoire*(real(3*ntheta+2,dp)*a2-a1)

cs = 1.0_dp-1.0_dp/(2.0_dp*aMoire)
sn = sqrt(1.0_dp-cs**2)
RotMatrix = reshape([cos(0.5_dp*acos(cs)),-sin(0.5_dp*acos(cs)),sin(0.5_dp*acos(cs)),cos(0.5_dp*acos(cs))],[2,2])

g1  = matmul(RotMatrix,g1)
g12 = matmul(RotMatrix,g12)

! fine grid of the band structure, independent of the numk of the Fock matrix

! Moire lattice parameters

t1 = real(ntheta,dp)*a1 + real(ntheta+1,dp)*a2
t2 = real(-ntheta-1,dp)*a1 + real(2*ntheta+1,dp)*a2
t3 = t2 - t1

call WignerSeitzCell(Coords,t1,t2,cs,sn)

do n=1,ndim
    Coords(n,1:2) = matmul(RotMatrix,Coords(n,1:2))
end do

t1 = matmul(RotMatrix,t1)
t2 = matmul(RotMatrix,t2)
t3 = t2 - t1

if(nrelax.EQ.1)then
    if(nlayers.EQ.2)then
        call LatticeRelaxationKoshino(Coords,g1,g12)
    endif
endif
if(nrelax.EQ.2)then
    if(nlayers.EQ.2)then
        call LatticeRelaxationCarr(Coords,g1,g12)
    endif
endif
