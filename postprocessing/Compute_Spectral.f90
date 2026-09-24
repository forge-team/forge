! Cleaned from plotSpectral.f90; original modules retained.
! Input zFock is the one-particle reduced density matrix from FORGE.
! Geometry/input: FORGE dev; response kernels: Hamiltonian.f90. See README.md.

program compute_spectral

use omp_lib
use Hamiltonian
use PostProcessingInput, only: BuildGeometry, ReadFock, InputParameters, ValidateResponseModel
use Geometry, only: NumberNeighborCells, OrderNeighborCells
use Setup, only: alphaInput => alpha, alphaHInput => alphaH
use lapack_routines

 implicit none

 character(1024) :: filename
 character(220) :: parameters
 integer(dp) :: ivk1, ivk2, n, m, i, j, imu

 integer(dp) :: icount, ncount, ncount_dB, numkPlot, nbounds(4)
 integer(dp) :: ndimEV, nNP, il, iu, nband, nspin, it
 integer(dp) :: ind
 integer(dp) :: nOverlap
 real(dp) :: fmu, energy, fmuX, fmuY, fsum
 real(dp) :: efermi, CBand
 real(dp) :: DiamagX(nmu,nT), DiamagY(nmu,nT), fKinetic(nmu,nT), MagMomentX(nmu,nT), MagMomentY(nmu,nT), AreaMoire

 real(dp) :: vk(2), t1(2), t2(2), t3(2), vq1(2), vq12(2)
 real(dp) :: ang_mat(2,2), eDirac(1)

 real(dp) :: om, T
 integer(dp), allocatable :: npointsBZ(:,:), npointsBack(:,:)
 integer(dp), allocatable :: npointsBZ_dB(:,:), npointsBack_dB(:,:)
 integer(dp), allocatable :: ind_j(:), ind_l(:)
 real(dp), allocatable :: coord(:,:), energia(:,:), vkBZ(:,:,:)
 complex(dp), allocatable :: zh(:,:,:), zevTB(:,:), zvxy1(:,:), zvxy2(:,:), zvxyI(:,:), zvxyD1(:,:), zvxyD2(:,:)
 complex(dp), allocatable :: zvxyDI(:,:)

 real(dp), allocatable :: dos(:), temp1(:), temp2(:), temp3(:)
 real(dp), allocatable :: drudeTotX(:), drudeMagX(:), drudeChiX(:)
 real(dp), allocatable :: drudeTotY(:), drudeMagY(:), drudeChiY(:)

 real(dp), allocatable :: dens(:,:)
 real(dp), allocatable :: potential(:,:)
 real(dp), allocatable :: fLongRange(:,:)
 complex(dp), allocatable :: zfock(:,:,:,:)

 real(dp) :: OmegaMin, OmegaMax
 integer(dp), parameter :: nspacing=10000
 integer(dp) :: nomega, nshift
 real(dp), parameter :: Ac=4._dp*pi*pi ! normalization

 real(dp), allocatable :: vEnergy(:,:,:), vMomentum(:,:,:)
 real(dp), allocatable :: vVx1(:,:,:), vVy1(:,:,:), vVDx(:,:,:), vVDy(:,:,:)
 complex(dp), allocatable :: zVmx1(:,:,:,:), zVmy1(:,:,:,:)
 real(dp), allocatable :: vVx2(:,:,:), vVy2(:,:,:)
 complex(dp), allocatable :: zVmx2(:,:,:,:), zVmy2(:,:,:,:)
 real(dp), allocatable :: vVxI(:,:,:), vVyI(:,:,:)
 complex(dp), allocatable :: zVmxI(:,:,:,:), zVmyI(:,:,:,:)

 real(dp), allocatable :: sigmaTotX(:,:,:), sigmaMagX(:,:,:), sigmaChiX(:,:,:)
 real(dp), allocatable :: sigmaTotY(:,:,:), sigmaMagY(:,:,:), sigmaChiY(:,:,:)
 real(dp), allocatable :: pivotEnergy(:,:), pivotMomentum(:,:)
 real(dp), allocatable :: pivotVx1(:,:), pivotVy1(:,:)
 complex(dp), allocatable :: pivotVmx1(:,:,:), pivotVmy1(:,:,:)
 real(dp), allocatable :: pivotVx2(:,:), pivotVy2(:,:)
 complex(dp), allocatable :: pivotVmx2(:,:,:), pivotVmy2(:,:,:)
 real(dp), allocatable :: pivotVxI(:,:), pivotVyI(:,:)
 complex(dp), allocatable :: pivotVmxI(:,:,:), pivotVmyI(:,:,:)

 real(dp) :: UPlot

 real(dp) :: alpha, alphaH

 call ValidateResponseModel()

 ndimEV=numb
il=(ndim-2*ndimEV+2*ncb)/2+1
iu=(ndim+2*ncb)/2
nNP=(ndim/2-il)+1

numkPlot=6

ncount_dB=(numkPlot+1)*(numkPlot+1)

allocate(vkBZ(1:numkPlot+1,1:numkPlot+1,1:2))

allocate(dens(1:ndim,1:numS))
allocate(energia(ndimEV,1:numS))
allocate(coord(ndim,3))
allocate(zh(ndim,ndim,1:numS))
allocate(zevTB(ndim,ndim))
allocate(zvxy1(ndim,ndim),zvxy2(ndim,ndim),zvxyI(ndim,ndim))
allocate(zvxyD1(ndim,ndim),zvxyD2(ndim,ndim),zvxyDI(ndim,ndim))

allocate(npointsBZ_dB(ncount_dB,2))
allocate(npointsBack_dB(numkPlot+1,numkPlot+1))

 allocate(zVmx1(ndimEV,ndimEV,numkPlot+1,numkPlot+1),zVmy1(ndimEV,ndimEV,numkPlot+1,numkPlot+1))
 allocate(zVmx2(ndimEV,ndimEV,numkPlot+1,numkPlot+1),zVmy2(ndimEV,ndimEV,numkPlot+1,numkPlot+1))
 allocate(zVmxI(ndimEV,ndimEV,numkPlot+1,numkPlot+1),zVmyI(ndimEV,ndimEV,numkPlot+1,numkPlot+1))

 allocate(vEnergy(ndimEV,numkPlot+1,numkPlot+1))
 allocate(vMomentum(2,numkPlot+1,numkPlot+1))

 allocate(pivotEnergy(ndimEV,4),pivotMomentum(2,4))
 allocate(pivotVmx1(ndimEV,ndimEV,4),pivotVmy1(ndimEV,ndimEV,4))
 allocate(pivotVmx2(ndimEV,ndimEV,4),pivotVmy2(ndimEV,ndimEV,4))
 allocate(pivotVmxI(ndimEV,ndimEV,4),pivotVmyI(ndimEV,ndimEV,4))

 allocate(vVDx(ndimEV,numkPlot+1,numkPlot+1))
 allocate(vVDy(ndimEV,numkPlot+1,numkPlot+1))

allocate(potential(1:ndim,1:numS))
allocate(fLongRange(1:ndim,1:ndim))

 write(*,*) tperp,dp,numI

! Retained cells: same ordering as Main_BandStructure.
call NumberNeighborCells(numI, ind)
ncount = numkPlot*numkPlot
allocate(zFock(ndim,ndim,ind,numS))
zFock = cmplx(0.0_dp,0.0_dp,dp)
allocate(ind_j(ind),ind_l(ind))
call OrderNeighborCells(numI, ind, ind_j, ind_l)

! Geometry: shared implementation of Main_BandStructure's construction.
call BuildGeometry(coord, t1, t2, t3, vq1, vq12, ang_mat)

! The long-range interaction uses the shared dev screening routine.
call longRange(fLongRange,coord,ndim,t1,t2)

UPlot=U

! Get Points

nOverlap=4

nbounds(1)=1 ! min for i_v
nbounds(2)=ndimEV-ncb+nOverlap ! max for i_v
nbounds(3)=ndimEV-ncb+1-nOverlap ! min for i_c
nbounds(4)=ndimEV ! max for i_c

write(*,*) 'Test',dir

! Input state: filename and packed complex ordering from FORGE dev.
alpha = alphaInput
alphaH = alphaHInput
call InputParameters(parameters)
call ReadFock(zFock, ind, parameters)
do nspin=1,numS
    do n=1,ndim
        dens(n,nspin) = real(zFock(n,n,1,nspin),dp)-0.5_dp
    enddo
enddo

write(*,*) 'ReadFock-finished'

potential=0.0_dp

if(numS.EQ.1)then
do n=1,ndim
do m=1,ndim
potential(n,1)=potential(n,1)+2._dp*dens(m,1)*alphaH*fLongRange(n,m)
enddo
enddo
potential(:,1)=potential(:,1)+UPlot*dens(:,1)
endif

if(numS.EQ.2)then
do n=1,ndim
do m=1,ndim
potential(n,1)=potential(n,1)+(dens(m,1)+dens(m,2))*alphaH*fLongRange(n,m)
enddo
enddo

potential(:,2)=potential(:,1)+UPlot*dens(:,1)
potential(:,1)=potential(:,1)+UPlot*dens(:,2)
endif

! Start Initialization

write(*,*) 'Electronic Density of one band in 10^12 cm^-2',4._dp/(sqrt(3._dp)/2._dp*ndim/4._dp*2.46*2.46)*10000._dp ! unit of density of one

write(*,*) 'Parameters'

write(*,*) 'Precision',dp

write(*,*) 'Reduced Dimension',ndimEV

write(*,*) 'Active Bands',il,iu

write(*,*) 'Number of k-points',ncount

write(*,*) 'Get Parameters for Output'

ncount=numkPlot*numkPlot
allocate(npointsBZ(ncount,8))
allocate(npointsBack(numkPlot,numkPlot))

call samplePoints1(npointsBZ,npointsBack,numkPlot,ncount)

call samplePoints_dB(npointsBZ_dB,npointsBack_dB,numkPlot,ncount_dB)

do ivk1=0,numkPlot
do ivk2=0,numkPlot

vkBZ(ivk1+1,ivk2+1,1)=(ivk1*vq1(1)+ivk2*vq12(1))/real(numkPlot,dp)
vkBZ(ivk1+1,ivk2+1,2)=(ivk1*vq1(2)+ivk2*vq12(2))/real(numkPlot,dp)

enddo
enddo

nspin=1
write(*,*) 'Start Loop'

!$omp parallel do &
!$omp private(vk,icount,zh,energia,zvxy1,zvxy2,zvxyI,zevTB,zvxyD1,zvxyD2,zvxyDI) &
!$omp shared(t1,t2,t3,il,iu,ncount_dB,coord,alpha,potential,zfock,ind_j,ind_l,npointsBZ_dB) &
!$omp shared(vEnergy,vMomentum,zVmx1,zVmy1,zVmx2,zVmy2,zVmxI,zVmyI,vVDx,vVDy)
 do icount = 1,ncount_dB

    vk(:)=vkBZ(npointsBZ_dB(icount,1)+1,npointsBZ_dB(icount,2)+1,:)

   call ComplexHk_MF(zh(:,:,nspin),coord,potential(:,nspin),alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]))

   call diagonalize(zh(:,:,nspin),energia(:,nspin),'V',il,iu)

   vEnergy(:,npointsBZ_dB(icount,1)+1,npointsBZ_dB(icount,2)+1)=energia(:,nspin)
   vMomentum(:,npointsBZ_dB(icount,1)+1,npointsBZ_dB(icount,2)+1)=vk(:)

   zevTB(:,il:iu)=zh(:,1:ndimEV,nspin)

   call vxy12I_MF(zvxy1,zvxy2,zvxyI,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),1_dp)

   call matrix_elements(zvxy1,zevTB,zVmx1(:,:,npointsBZ_dB(icount,1)+1,npointsBZ_dB(icount,2)+1),il,iu,il,iu)

   call matrix_elements(zvxy2,zevTB,zVmx2(:,:,npointsBZ_dB(icount,1)+1,npointsBZ_dB(icount,2)+1),il,iu,il,iu)

   call matrix_elements(zvxyI,zevTB,zVmxI(:,:,npointsBZ_dB(icount,1)+1,npointsBZ_dB(icount,2)+1),il,iu,il,iu)

   call vxy12I_MF(zvxy1,zvxy2,zvxyI,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),2_dp)

   call matrix_elements(zvxy1,zevTB,zVmy1(:,:,npointsBZ_dB(icount,1)+1,npointsBZ_dB(icount,2)+1),il,iu,il,iu)

   call matrix_elements(zvxy2,zevTB,zVmy2(:,:,npointsBZ_dB(icount,1)+1,npointsBZ_dB(icount,2)+1),il,iu,il,iu)

   call matrix_elements(zvxyI,zevTB,zVmyI(:,:,npointsBZ_dB(icount,1)+1,npointsBZ_dB(icount,2)+1),il,iu,il,iu)

   call vxy12ID_MF(zvxyD1,zvxyD2,zvxyDI,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),1_dp)

   call z_matrix_elements_diag_real(zvxyD1+zvxyD2+zvxyDI,zevTB,vVDx(:,npointsBZ_dB(icount,1)+1,npointsBZ_dB(icount,2)+1))

   call vxy12ID_MF(zvxyD1,zvxyD2,zvxyDI,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),2_dp)

   call z_matrix_elements_diag_real(zvxyD1+zvxyD2+zvxyDI,zevTB,vVDy(:,npointsBZ_dB(icount,1)+1,npointsBZ_dB(icount,2)+1))

enddo
!$omp end parallel do

write(*,*) 'End Loop'

deallocate(zh)
deallocate(zvxy1,zvxy2,zvxyI)

deallocate(potential)
deallocate(fLongRange)

! Set spectral range

efermi=(vEnergy(nNP,numkPlot/3+1,numkPlot/3+1)+vEnergy(nNP+1,numkPlot/3+1,numkPlot/3+1))/2._dp

write(*,*) efermi

fmu=efermi
eDirac(nspin)=fmu

write(*,*) vEnergy(nNP,numkPlot/3+1,numkPlot/3+1),vEnergy(nNP+1,numkPlot/3+1,numkPlot/3+1)

write(*,*) fmu,nNP,numkPlot

write(*,*) maxval(vEnergy(ndimEV,:,:)),minval(vEnergy(1,:,:))

 CBand=maxval(vEnergy(nNP+1,:,:))-minval(vEnergy(nNP+1,:,:))

 CBand=CBand*5

  OmegaMax=maxval(vEnergy(ndimEV,:,:))-fmu
  OmegaMin=minval(vEnergy(1,:,:))-fmu

  OmegaMax=OmegaMax*1.1_dp
  OmegaMin=OmegaMin*1.1_dp

  write(*,*) 'maxmin',OmegaMax,OmegaMin

  nomega=int((OmegaMax-OmegaMin)*nspacing)
  nshift=int(-OmegaMin*nspacing)

  write(*,*) OmegaMax,OmegaMin
  write(*,*) nomega,nshift

! Diamagnetic Current

 DiamagX=0._dp
 DiamagY=0._dp
 fKinetic=0._dp
 MagMomentX=0._dp
 MagMomentY=0._dp

AreaMoire=sqrt(3._dp)/2._dp*real(numkPlot*numkPlot,dp)*(3._dp*ntheta*ntheta+3._dp*ntheta+1._dp)

do imu=1,nmu

fmu=eDirac(nspin)!    CBand*(nmu-imu)/real(nmu,8)

do iT=1,nT ! iT

T=ftemperature(iT)!*0.0000861733_dp ! Temperature

fsum=0._dp

   do icount = 1,ncount

    do nband=1,ndimEV

    energy=vEnergy(nband,npointsBZ(icount,1)+1,npointsBZ(icount,2)+1)-fmu

    fmuX=real(zVmy1(nband,nband,npointsBZ(icount,1)+1,npointsBZ(icount,2)+1))-real(zVmy2(nband,nband,npointsBZ(icount,1)+1,npointsBZ(icount,2)+1))
    fmuY=real(zVmx2(nband,nband,npointsBZ(icount,1)+1,npointsBZ(icount,2)+1))-real(zVmx1(nband,nband,npointsBZ(icount,1)+1,npointsBZ(icount,2)+1))

    DiamagX(imu,iT)=DiamagX(imu,iT)+vVDx(nband,npointsBZ(icount,1)+1,npointsBZ(icount,2)+1)*fermi_dist(energy,T)
    DiamagY(imu,iT)=DiamagY(imu,iT)+vVDy(nband,npointsBZ(icount,1)+1,npointsBZ(icount,2)+1)*fermi_dist(energy,T)
    fKinetic(imu,iT)=fKinetic(imu,iT)+vEnergy(nband,npointsBZ(icount,1)+1,npointsBZ(icount,2)+1)*fermi_dist(energy,T)
    MagMomentX(imu,iT)=MagMomentX(imu,iT)+fmuX*fermi_dist(energy,T)
    MagMomentY(imu,iT)=MagMomentY(imu,iT)+fmuY*fermi_dist(energy,T)

    enddo

    do nband=1,ndimEV/2_dp

    fsum=fsum+vEnergy(nband,npointsBZ(icount,1)+1,npointsBZ(icount,2)+1)

    enddo

    enddo

    write(*,*) imu,iT,DiamagX(imu,iT),DiamagX(imu,iT)/AreaMoire
    write(*,*) imu,iT,DiamagY(imu,iT),DiamagY(imu,iT)/AreaMoire
    write(*,*) imu,iT,fKinetic(imu,iT)/6._dp,fsum/6._dp,fKinetic(imu,iT)/AreaMoire/6._dp
    write(*,*) 'Magnetic Moment'
    write(*,*) imu,iT,MagMomentX(imu,iT)*137._dp/300._dp*3.4/0.529/sqrt(3.)*2._dp/AreaMoire
    write(*,*) imu,iT,MagMomentY(imu,iT)*137._dp/300._dp*3.4/0.529/sqrt(3.)*2._dp/AreaMoire

 enddo
 enddo

    DiamagX=DiamagX/AreaMoire
    DiamagY=DiamagY/AreaMoire

 write(*,*) 'Test Output'

write(filename,'(A,A5,A)') trim(dir),'Bands',trim(parameters)
open(99,file=filename,status='replace')

call plotBandsK(vEnergy,fmu,npointsBack,ndimEV,numkPlot,ncount,vq1,vq12)

write(filename,'(A,A7,A)') trim(dir),'BandsDP',trim(parameters)
open(99,file=filename,status='replace')

nband=int(nfilling/2.)

call plotBandsDensity_dB(vEnergy,npointsBZ_dB,ndimEV,numkPlot,nband,ncount_dB,nNP)

write(*,*) 'Test Drude'

if(nDrude.EQ.1.OR.nDos.EQ.1)then

allocate(vVx1(ndimEV,numkPlot+1,numkPlot+1))
allocate(vVy1(ndimEV,numkPlot+1,numkPlot+1))
allocate(pivotVx1(ndimEV,4),pivotVy1(ndimEV,4))
allocate(vVx2(ndimEV,numkPlot+1,numkPlot+1))
allocate(vVy2(ndimEV,numkPlot+1,numkPlot+1))
allocate(pivotVx2(ndimEV,4),pivotVy2(ndimEV,4))
allocate(vVxI(ndimEV,numkPlot+1,numkPlot+1))
allocate(vVyI(ndimEV,numkPlot+1,numkPlot+1))
allocate(pivotVxI(ndimEV,4),pivotVyI(ndimEV,4))
allocate(temp1(nomega))
allocate(temp2(nomega))
allocate(temp3(nomega))
allocate(drudeTotX(nomega))
allocate(drudeMagX(nomega))
allocate(drudeChiX(nomega))
allocate(drudeTotY(nomega))
allocate(drudeMagY(nomega))
allocate(drudeChiY(nomega))
allocate(dos(nomega))

do icount=1,ncount_dB
ivk1=npointsBZ_dB(icount,1)
ivk2=npointsBZ_dB(icount,2)

do nband=1,ndimEV
vVx1(nband,ivk1+1,ivk2+1)=real(zVmx1(nband,nband,ivk1+1,ivk2+1))
vVy1(nband,ivk1+1,ivk2+1)=real(zVmy1(nband,nband,ivk1+1,ivk2+1))
vVx2(nband,ivk1+1,ivk2+1)=real(zVmx2(nband,nband,ivk1+1,ivk2+1))
vVy2(nband,ivk1+1,ivk2+1)=real(zVmy2(nband,nband,ivk1+1,ivk2+1))
vVxI(nband,ivk1+1,ivk2+1)=real(zVmxI(nband,nband,ivk1+1,ivk2+1))
vVyI(nband,ivk1+1,ivk2+1)=real(zVmyI(nband,nband,ivk1+1,ivk2+1))
enddo

enddo

endif

if(nDrude.EQ.1)then

write(*,*) 'Test Drude2'

drudeTotX=0._dp
drudeMagX=0._dp
drudeChiX=0._dp
drudeTotY=0._dp
drudeMagY=0._dp
drudeChiY=0._dp

temp1=0._dp
temp2=0._dp
temp3=0._dp

!$omp parallel do &
!$omp private(icount,ivk1,ivk2,pivotEnergy,pivotMomentum,pivotVx1,pivotVy1,pivotVx2,pivotVy2,pivotVxI,pivotVyI) &
!$omp private(temp1,temp2,temp3) &
!$omp shared(ncount,vEnergy,vMomentum,npointsBZ,vVx1,vVy1,vVx2,vVy2,vVxI,vVyI) &
!$omp shared(nomega,nshift,fmu,nbounds) &
!$omp reduction(+:drudeTotX) &
!$omp reduction(+:drudeMagX) &
!$omp reduction(+:drudeChiX) &
!$omp reduction(+:drudeTotY) &
!$omp reduction(+:drudeMagY) &
!$omp reduction(+:drudeChiY)

do icount=1,ncount

ivk1=npointsBZ(icount,1)
ivk2=npointsBZ(icount,2)

pivotEnergy=reshape([vEnergy(:,ivk1+1,ivk2+1),vEnergy(:,ivk1+1,ivk2+2),vEnergy(:,ivk1+2,ivk2+1),vEnergy(:,ivk1+2,ivk2+2)],[ndimEV,4_dp])

pivotVx1=reshape([vVx1(:,ivk1+1,ivk2+1),vVx1(:,ivk1+1,ivk2+2),vVx1(:,ivk1+2,ivk2+1),vVx1(:,ivk1+2,ivk2+2)],[ndimEV,4_dp])

pivotVy1=reshape([vVy1(:,ivk1+1,ivk2+1),vVy1(:,ivk1+1,ivk2+2),vVy1(:,ivk1+2,ivk2+1),vVy1(:,ivk1+2,ivk2+2)],[ndimEV,4_dp])

pivotVx2=reshape([vVx2(:,ivk1+1,ivk2+1),vVx2(:,ivk1+1,ivk2+2),vVx2(:,ivk1+2,ivk2+1),vVx2(:,ivk1+2,ivk2+2)],[ndimEV,4_dp])

pivotVy2=reshape([vVy2(:,ivk1+1,ivk2+1),vVy2(:,ivk1+1,ivk2+2),vVy2(:,ivk1+2,ivk2+1),vVy2(:,ivk1+2,ivk2+2)],[ndimEV,4_dp])

pivotVxI=reshape([vVxI(:,ivk1+1,ivk2+1),vVxI(:,ivk1+1,ivk2+2),vVxI(:,ivk1+2,ivk2+1),vVxI(:,ivk1+2,ivk2+2)],[ndimEV,4_dp])

pivotVyI=reshape([vVyI(:,ivk1+1,ivk2+1),vVyI(:,ivk1+1,ivk2+2),vVyI(:,ivk1+2,ivk2+1),vVyI(:,ivk1+2,ivk2+2)],[ndimEV,4_dp])

pivotMomentum=reshape([vMomentum(:,ivk1+1,ivk2+1),vMomentum(:,ivk1+1,ivk2+2),vMomentum(:,ivk1+2,ivk2+1),vMomentum(:,ivk1+2,ivk2+2)],[2_dp,4_dp])

call drudeFullX(temp1,temp2,temp3,pivotMomentum,pivotEnergy,pivotVx1,pivotVy1,pivotVx2,pivotVy2,pivotVxI,pivotVyI,nomega,nspacing,nshift,fmu,nbounds)

drudeTotX(:)=drudeTotX(:)+temp1(:)
drudeMagX(:)=drudeMagX(:)+temp2(:)
drudeChiX(:)=drudeChiX(:)+temp3(:)

call drudeFullY(temp1,temp2,temp3,pivotMomentum,pivotEnergy,pivotVx1,pivotVy1,pivotVx2,pivotVy2,pivotVxI,pivotVyI,nomega,nspacing,nshift,fmu,nbounds)

drudeTotY(:)=drudeTotY(:)+temp1(:)
drudeMagY(:)=drudeMagY(:)+temp2(:)
drudeChiY(:)=drudeChiY(:)+temp3(:)

enddo ! icount
!$omp end parallel do

write(filename,'(A,A17,I0,A9,I0,A)') trim(dir),'DrudeTot-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(87,file=filename,status='replace')

write(filename,'(A,A17,I0,A9,I0,A)') trim(dir),'DrudeMag-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(88,file=filename,status='replace')

write(filename,'(A,A17,I0,A9,I0,A)') trim(dir),'DrudeChi-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(89,file=filename,status='replace')

write(filename,'(A,A18,I0,A9,I0,A)') trim(dir),'DrudeTotD-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(97,file=filename,status='replace')

write(filename,'(A,A18,I0,A9,I0,A)') trim(dir),'DrudeMagD-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(98,file=filename,status='replace')

write(filename,'(A,A18,I0,A9,I0,A)') trim(dir),'DrudeChiD-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(99,file=filename,status='replace')

! Normalization to "sigma_0=e^2/hbar/8" - spin not included
! Implementation with "exact" symmetry for negative frequencies

do n=1,nomega
om=real(n-nshift)/real(nspacing)
write(87,*) om,drudeTotX(n)/2._dp/Ac+drudeTotY(n)/2._dp/Ac !spin is not included in units of (e/hbar)^2
write(88,*) om,drudeMagX(n)/2._dp/Ac+drudeMagY(n)/2._dp/Ac !spin is not included in units of (e/hbar)^2
write(89,*) om,drudeChiX(n)/2._dp/Ac+drudeChiY(n)/2._dp/Ac !spin is not included in units of (e/hbar)^2
write(97,*) om,drudeTotX(n)/Ac,drudeTotY(n)/Ac !spin is not included in units of (e/hbar)^2
write(98,*) om,drudeMagX(n)/Ac,drudeMagY(n)/Ac !spin is not included in units of (e/hbar)^2
write(99,*) om,drudeChiX(n)/Ac,drudeChiY(n)/Ac !spin is not included in units of (e/hbar)^2
enddo
close(87)
close(88)
close(89)
close(97)
close(98)
close(99)

endif

write(*,*) 'DOS'

if(nDos.EQ.1)then

dos=0._dp
temp1=0._dp

!$omp parallel do &
!$omp private(icount,ivk1,ivk2,pivotEnergy,pivotMomentum) &
!$omp private(temp1) &
!$omp shared(ncount,vEnergy,vMomentum,npointsBZ,nomega,nshift,fmu,nbounds) &
!$omp reduction(+:dos)
do icount=1,ncount

ivk1=npointsBZ(icount,1)
ivk2=npointsBZ(icount,2)

pivotEnergy=reshape([vEnergy(:,ivk1+1,ivk2+1),vEnergy(:,ivk1+1,ivk2+2),vEnergy(:,ivk1+2,ivk2+1),vEnergy(:,ivk1+2,ivk2+2)],[ndimEV,4_dp])

pivotMomentum=reshape([vMomentum(:,ivk1+1,ivk2+1),vMomentum(:,ivk1+1,ivk2+2),vMomentum(:,ivk1+2,ivk2+1),vMomentum(:,ivk1+2,ivk2+2)],[2_dp,4_dp])

call dos_F(temp1,pivotMomentum,pivotEnergy,nomega,nspacing,nshift,fmu,nbounds)

dos(:)=dos(:)+temp1(:)

enddo ! icount
!$omp end parallel do

write(*,*) 'WriteDos'

write(filename,'(A,A12,I0,A9,I0,A)') trim(dir),'Dos-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(99,file=filename,status='replace')

write(filename,'(A,A18,I0,A9,I0,A)') trim(dir),'DrudeMaki-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(79,file=filename,status='replace')

write(filename,'(A,A21,I0,A9,I0,A)') trim(dir),'DrudeMakiChi-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(89,file=filename,status='replace')

do n=1,nomega
om=real(n-nshift)/real(nspacing)
write(99,*) om,dos(n)/Ac,2.*abs(om)/pi/(3./4.)/2.7/2.7 ! Spin is not included
if(dos(n).GT.0.)then
write(79,*) om,(drudeMagX(n)/2._dp+drudeMagY(n)/2._dp)/dos(n) !spin is not included in units of (e/hbar)^2
else
write(79,*) om,0._dp
endif
if(dos(n).GT.0.)then
write(89,*) om,(drudeChiX(n)/2._dp+drudeChiY(n)/2._dp)/dos(n) !spin is not included in units of (e/hbar)^2
else
write(89,*) om,0._dp
endif
enddo

! normalization: int DOS*sqrt(3)/2*ndim/4= number of bands =iu-il

close(99)
close(79)
close(89)

deallocate(dos)
deallocate(temp1)
deallocate(temp2)
deallocate(temp3)

endif

  OmegaMax=maxval(vEnergy(ndimEV,:,:))-minval(vEnergy(1,:,:))

  OmegaMin=-OmegaMax*0.1_dp
  OmegaMax=OmegaMax*1.1_dp

  write(*,*) 'maxmin',OmegaMax,OmegaMin

  nomega=int((OmegaMax-OmegaMin)*nspacing)
  nshift=int(-OmegaMin*nspacing)

  write(*,*) OmegaMax,OmegaMin
  write(*,*) nomega,nshift

  allocate(sigmaTotX(nomega,nmu,nT))
  allocate(sigmaMagX(nomega,nmu,nT))
  allocate(sigmaChiX(nomega,nmu,nT))
  allocate(sigmaTotY(nomega,nmu,nT))
  allocate(sigmaMagY(nomega,nmu,nT))
  allocate(sigmaChiY(nomega,nmu,nT))
  allocate(temp1(nomega))
  allocate(temp2(nomega))
  allocate(temp3(nomega))

sigmaTotX=0._dp
sigmaMagX=0._dp
sigmaChiX=0._dp
sigmaTotY=0._dp
sigmaMagY=0._dp
sigmaChiY=0._dp

 if(nfilling.EQ.0)then
 CBand=maxval(vEnergy(nNP+2,:,:))-minval(vEnergy(nNP+1,:,:))
 endif

 if(nfilling.EQ.2)then
 CBand=maxval(vEnergy(nNP+2,:,:))-minval(vEnergy(nNP+2,:,:))
 endif

 if(nfilling.EQ.-2)then
 CBand=maxval(vEnergy(nNP+2,:,:))-minval(vEnergy(nNP,:,:))
 endif

do imu=1,nmu

fmu=eDirac(nspin) + CBand*(nmu-imu)/real(nmu,8)

do iT=1,nT ! iT

T=ftemperature(iT)!*0.0000861733_dp ! Temperature

!$omp parallel do &
!$omp private(icount,ivk1,ivk2,pivotEnergy,pivotMomentum,pivotVmx1,pivotVmy1,pivotVmx2,pivotVmy2,pivotVmxI,pivotVmyI) &
!$omp private(temp1,temp2,temp3) &
!$omp shared(ncount,iT,fmu,T,vEnergy,vMomentum,zVmx1,zVmy1,zVmx2,zVmy2,zVmxI,zVmyI,nomega,nshift,nbounds) &
!$omp reduction(+:sigmaTotX) &
!$omp reduction(+:sigmaMagX) &
!$omp reduction(+:sigmaChiX) &
!$omp reduction(+:sigmaTotY) &
!$omp reduction(+:sigmaMagY) &
!$omp reduction(+:sigmaChiY)
do icount=1,ncount

ivk1=npointsBZ(icount,1)
ivk2=npointsBZ(icount,2)

pivotEnergy=reshape([vEnergy(:,ivk1+1,ivk2+1),vEnergy(:,ivk1+1,ivk2+2),vEnergy(:,ivk1+2,ivk2+1),vEnergy(:,ivk1+2,ivk2+2)],[ndimEV,4_dp])

pivotVmx1=reshape([zVmx1(:,:,ivk1+1,ivk2+1),zVmx1(:,:,ivk1+1,ivk2+2),zVmx1(:,:,ivk1+2,ivk2+1),zVmx1(:,:,ivk1+2,ivk2+2)],[ndimEV,ndimEV,4_dp])

pivotVmy1=reshape([zVmy1(:,:,ivk1+1,ivk2+1),zVmy1(:,:,ivk1+1,ivk2+2),zVmy1(:,:,ivk1+2,ivk2+1),zVmy1(:,:,ivk1+2,ivk2+2)],[ndimEV,ndimEV,4_dp])

pivotVmx2=reshape([zVmx2(:,:,ivk1+1,ivk2+1),zVmx2(:,:,ivk1+1,ivk2+2),zVmx2(:,:,ivk1+2,ivk2+1),zVmx2(:,:,ivk1+2,ivk2+2)],[ndimEV,ndimEV,4_dp])

pivotVmy2=reshape([zVmy2(:,:,ivk1+1,ivk2+1),zVmy2(:,:,ivk1+1,ivk2+2),zVmy2(:,:,ivk1+2,ivk2+1),zVmy2(:,:,ivk1+2,ivk2+2)],[ndimEV,ndimEV,4_dp])

pivotVmxI=reshape([zVmxI(:,:,ivk1+1,ivk2+1),zVmxI(:,:,ivk1+1,ivk2+2),zVmxI(:,:,ivk1+2,ivk2+1),zVmxI(:,:,ivk1+2,ivk2+2)],[ndimEV,ndimEV,4_dp])

pivotVmyI=reshape([zVmyI(:,:,ivk1+1,ivk2+1),zVmyI(:,:,ivk1+1,ivk2+2),zVmyI(:,:,ivk1+2,ivk2+1),zVmyI(:,:,ivk1+2,ivk2+2)],[ndimEV,ndimEV,4_dp])

pivotMomentum=reshape([vMomentum(:,ivk1+1,ivk2+1),vMomentum(:,ivk1+1,ivk2+2),vMomentum(:,ivk1+2,ivk2+1),vMomentum(:,ivk1+2,ivk2+2)],[2_dp,4_dp])

call conductivityFullX(temp1,temp2,temp3,pivotMomentum,pivotEnergy,pivotVmx1,pivotVmy1,pivotVmx2,pivotVmy2,pivotVmxI,pivotVmyI,nomega,nspacing,nshift,fmu,T,nbounds)

sigmaTotX(:,imu,iT)=sigmaTotX(:,imu,iT)+temp1(:)
sigmaMagX(:,imu,iT)=sigmaMagX(:,imu,iT)+temp2(:)
sigmaChiX(:,imu,iT)=sigmaChiX(:,imu,iT)+temp3(:)

call conductivityFullY(temp1,temp2,temp3,pivotMomentum,pivotEnergy,pivotVmx1,pivotVmy1,pivotVmx2,pivotVmy2,pivotVmxI,pivotVmyI,nomega,nspacing,nshift,fmu,T,nbounds)

sigmaTotY(:,imu,iT)=sigmaTotY(:,imu,iT)+temp1(:)
sigmaMagY(:,imu,iT)=sigmaMagY(:,imu,iT)+temp2(:)
sigmaChiY(:,imu,iT)=sigmaChiY(:,imu,iT)+temp3(:)

enddo ! icount
!$omp end parallel do

enddo ! iT

enddo ! imu

! Kramers-Kronig postprocessing: compute_kramers_kronig.f90.

write(*,*) 'Output'

do imu=1,nmu

fmu=CBand*(nmu-imu)/real(nmu,8)

fmu=0._dp!fmu=eDirac(nspin)

write(*,*) imu,fmu

do iT=1,nT ! iT

T=ftemperature(iT)*0.0000861733_dp

write(filename,'(A,A15,I0,A3,F5.3,A9,I0,A9,I0,A)') trim(dir),'ConductivityTot',int(ftemperature(iT)),'-Mu',fmu,'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(87,file=filename,status='replace')

write(filename,'(A,A15,I0,A3,F5.3,A9,I0,A9,I0,A)') trim(dir),'ConductivityMag',int(ftemperature(iT)),'-Mu',fmu,'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(88,file=filename,status='replace')

write(filename,'(A,A15,I0,A3,F5.3,A9,I0,A9,I0,A)') trim(dir),'ConductivityChi',int(ftemperature(iT)),'-Mu',fmu,'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(89,file=filename,status='replace')

write(filename,'(A,A16,I0,A3,F5.3,A9,I0,A9,I0,A)') trim(dir),'ConductivityTotD',int(ftemperature(iT)),'-Mu',fmu,'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(97,file=filename,status='replace')

write(filename,'(A,A16,I0,A3,F5.3,A9,I0,A9,I0,A)') trim(dir),'ConductivityMagD',int(ftemperature(iT)),'-Mu',fmu,'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(98,file=filename,status='replace')

write(filename,'(A,A16,I0,A3,F5.3,A9,I0,A9,I0,A)') trim(dir),'ConductivityChiD',int(ftemperature(iT)),'-Mu',fmu,'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(99,file=filename,status='replace')

! Normalization to "sigma_0=e^2/hbar/8"
! Implementation with "exact" symmetry for negative frequencies
do n=1,nomega
om=real(n-nshift)/real(nspacing)
write(87,*) om,sigmaTotX(n,imu,iT)*8_dp/Ac/2._dp+sigmaTotY(n,imu,iT)*8_dp/Ac/2._dp
write(88,*) om,sigmaMagX(n,imu,iT)*8_dp/Ac/2._dp+sigmaMagY(n,imu,iT)*8_dp/Ac/2._dp
write(89,*) om,sigmaChiX(n,imu,iT)*8_dp/Ac/2._dp+sigmaChiY(n,imu,iT)*8_dp/Ac/2._dp
write(97,*) om,sigmaTotX(n,imu,iT)*8_dp/Ac,sigmaTotY(n,imu,iT)*8_dp/Ac
write(98,*) om,sigmaMagX(n,imu,iT)*8_dp/Ac,sigmaMagY(n,imu,iT)*8_dp/Ac
write(99,*) om,sigmaChiX(n,imu,iT)*8_dp/Ac,sigmaChiY(n,imu,iT)*8_dp/Ac
enddo
close(87)
close(88)
close(89)
close(97)
close(98)
close(99)

enddo !iT

enddo !imu

write(*,*) 'Test'

write(filename,'(A,A5,A)') trim(dir),'CBand',trim(parameters)
open(87,file=filename,status='replace')

write(87,*) CBand,numkPlot,nspacing

close(87)

do imu=1,nmu

fmu=CBand*(nmu-imu)/real(nmu,8)

write(*,*) imu,fmu

do iT=1,nT ! iT

T=ftemperature(iT)*0.0000861733_dp

write(filename,'(A,A6,I0,A3,F5.3,A9,I0,A9,I0,A)') trim(dir),'DiaTot',int(ftemperature(iT)),'-Mu',fmu,'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(67,file=filename,status='replace')

write(filename,'(A,A7,I0,A3,F5.3,A9,I0,A9,I0,A)') trim(dir),'DiaTotX',int(ftemperature(iT)),'-Mu',fmu,'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(76,file=filename,status='replace')

write(filename,'(A,A7,I0,A3,F5.3,A9,I0,A9,I0,A)') trim(dir),'DiaTotY',int(ftemperature(iT)),'-Mu',fmu,'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(77,file=filename,status='replace')

! Normalization to "sigma_0=e^2/hbar/8"
! Implementation with "exact" symmetry for negative frequencies

write(67,*) nomega,nshift,(DiamagX(imu,iT)+DiamagY(imu,iT))/2._dp

write(76,*) nomega,nshift,DiamagX(imu,iT)

write(77,*) nomega,nshift,DiamagY(imu,iT)

close(67)
close(87)
close(77)
close(97)

enddo !iT

enddo !imu

write(*,*) 'TestEnd'

end program compute_spectral
