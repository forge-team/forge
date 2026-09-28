! Cleaned from plotSpectralOptimize.f90; original modules retained.
! Input zFock is the one-particle reduced density matrix from FORGE.
! Geometry/input: FORGE dev; response kernels: Hamiltonian.f90. See README.md.

program compute_spectral_optimize

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

 integer(dp) :: icount, ncount, numkPlot, nbounds(4)
 integer(dp) :: ndimEV, nNP, il, iu, nband, nspin, it
 integer(dp) :: ind
 integer(dp) :: nOverlap
 real(dp) :: fmu, sumSigma
 real(dp) :: efermi, CBand
 real(dp) :: DiamagX(nmu,nT), DiamagY(nmu,nT), fKinetic(nmu,nT), DiamagNO(nT)
 real(dp) :: AreaMoire

 real(dp) :: vk(2), t1(2), t2(2), t3(2), vq1(2), vq12(2)
 real(dp) :: ang_mat(2,2)

 real(dp) :: om, T
 integer(dp), allocatable :: ind_j(:), ind_l(:)
 real(dp), allocatable :: coord(:,:), energia(:,:), vkBZ(:,:,:)
 complex(dp), allocatable :: zh(:,:,:), zevTB(:,:), zvxy1(:,:), zvxy2(:,:), zvxyI(:,:)

 real(dp), allocatable :: dos(:), drude(:), temp1(:), tempD(:)

 real(dp), allocatable :: dens(:,:)
 real(dp), allocatable :: potential(:,:)
 real(dp), allocatable :: fLongRange(:,:)
 complex(dp), allocatable :: zfock(:,:,:,:)

 real(dp), parameter :: OmegaMinD = -2.3_dp
 real(dp), parameter :: OmegaMaxD = 2.3_dp
 integer(dp), parameter :: nspacingD=10000
 integer(dp), parameter :: nomegaD=int((OmegaMaxD-OmegaMinD)*nspacingD)
 integer(dp), parameter :: nshiftD=int(-OmegaMinD*nspacingD)

 real(dp) :: OmegaMin, OmegaMax
 integer(dp), parameter :: nspacing=10000
 integer(dp) :: nomega, nshift
 real(dp), parameter :: Ac=4._dp*pi*pi ! normalization

 real(dp), allocatable :: vEnergy(:,:,:), vMomentum(:,:,:)
 real(dp), allocatable :: vVx(:,:,:), vVy(:,:,:), vVDx(:), vVDy(:), vVDxNO(:), vVDyNO(:)
 complex(dp), allocatable :: zVmx(:,:,:,:), zVmy(:,:,:,:)

 real(dp), allocatable :: sigmaTot(:,:,:)
 real(dp), allocatable :: pivotEnergy(:,:), pivotMomentum(:,:)
 real(dp), allocatable :: pivotVx(:,:), pivotVy(:,:)
 complex(dp), allocatable :: pivotVmx(:,:,:), pivotVmy(:,:,:)

 real(dp) :: UPlot

 real(dp) :: alpha, alphaH

 call ValidateResponseModel()

 ndimEV=numb
il=(ndim-2*ndimEV+2*ncb)/2+1
iu=(ndim+2*ncb)/2
nNP=(ndim/2-il)+1

numkPlot=30

allocate(vkBZ(1:numkPlot+1,1:numkPlot+1,1:2))

allocate(dens(1:ndim,1:numS))
allocate(energia(ndimEV,1:numS))
allocate(coord(ndim,3))
allocate(zh(ndim,ndim,1:numS))
allocate(zevTB(ndim,ndim))
allocate(zvxy1(ndim,ndim),zvxy2(ndim,ndim),zvxyI(ndim,ndim))

 allocate(zVmx(ndimEV,ndimEV,numkPlot+1,3),zVmy(ndimEV,ndimEV,numkPlot+1,3))
 allocate(vVx(ndimEV,numkPlot+1,3),vVy(ndimEV,numkPlot+1,3))

 allocate(vEnergy(ndimEV,numkPlot+1,3))
 allocate(vMomentum(2,numkPlot+1,3))

 allocate(pivotEnergy(ndimEV,4),pivotMomentum(2,4))

 allocate(vVDx(ndimEV))
 allocate(vVDy(ndimEV))
 allocate(vVDxNO(ndimEV))
 allocate(vVDyNO(ndimEV))

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

do ivk1=0,numkPlot
do ivk2=0,numkPlot

vkBZ(ivk1+1,ivk2+1,:)=(ivk1*vq1(:)+ivk2*vq12(:))/real(numkPlot,dp)

enddo
enddo

 DiamagX=0._dp
 DiamagY=0._dp
 fKinetic=0._dp

nspin=1

   vk(:)=vkBZ(numkPlot/3+1,numkPlot/3+1,:)

   call ComplexHk_MF(zh(:,:,nspin),coord,potential(:,nspin),alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]))

   call diagonalize(zh(:,:,nspin),energia(:,nspin),'N',il,iu)

efermi=(energia(nNP,nspin)+energia(nNP+1,nspin))/2._dp

fmu=efermi

write(*,*) efermi
write(*,*) -energia(nNP,nspin)+energia(nNP+1,nspin)

write(*,*) fmu,nNP,numkPlot

 CBand=energia(nNP+3,nspin)-energia(nNP+1,nspin)

 CBand=CBand*2

 write(*,*) CBand,energia(nNP+3,nspin),energia(nNP+1,nspin)

  OmegaMax=energia(ndimEV,nspin)-fmu
  OmegaMin=energia(1,nspin)-fmu

  write(*,*) 'maxmin',OmegaMax,OmegaMin

  OmegaMax=OmegaMax-OmegaMin

  OmegaMin=-OmegaMax*0.1_dp
  OmegaMax=OmegaMax*1.1_dp

  write(*,*) 'maxmin',OmegaMax,OmegaMin

  nomega=int((OmegaMax-OmegaMin)*nspacing)
  nshift=int(-OmegaMin*nspacing)

  write(*,*) OmegaMax,OmegaMin
  write(*,*) nomega,nshift

  allocate(sigmaTot(nomega,nmu,nT))
  allocate(dos(nomegaD))
  allocate(drude(nomegaD))
  allocate(temp1(nomega))
  allocate(tempD(nomegaD))

sigmaTot=0._dp
dos=0._dp
drude=0._dp

 DiamagX=0._dp
 DiamagY=0._dp
 DiamagNO=0._dp

AreaMoire=sqrt(3._dp)/2._dp*real(numkPlot*numkPlot,dp)*(3._dp*ntheta*ntheta+3._dp*ntheta+1._dp)

nspin=1
write(*,*) 'Start Loop'

ivk1=0
!$omp parallel do &
!$omp private(vk,ivk2,zh,energia,zvxy1,zvxy2,zvxyI,zevTB,nband) &
!$omp shared(ivk1,t1,t2,t3,il,iu,coord,alpha,potential,zfock,ind_j,ind_l,numkPlot,vkBZ) &
!$omp shared(vEnergy,vMomentum,zVmx,zVmy)
 do ivk2 = 0,numkPlot

    vk(:)=vkBZ(ivk1+1,ivk2+1,:)

   call ComplexHk_MF(zh(:,:,nspin),coord,potential(:,nspin),alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]))

   call diagonalize(zh(:,:,nspin),energia(:,nspin),'V',il,iu)

   vEnergy(:,ivk2+1,1)=energia(:,nspin)
   vMomentum(:,ivk2+1,1)=vk(:)

   zevTB(:,il:iu)=zh(:,1:ndimEV,nspin)

   call vxy12I_MF(zvxy1,zvxy2,zvxyI,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),1_dp)

   call matrix_elements(zvxy1+zvxy2+zvxyI,zevTB,zVmx(:,:,ivk2+1,1),il,iu,il,iu)

   call vxy12I_MF(zvxy1,zvxy2,zvxyI,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),2_dp)

   call matrix_elements(zvxy1+zvxy2+zvxyI,zevTB,zVmy(:,:,ivk2+1,1),il,iu,il,iu)

   vVx(:,ivk2+1,1)=[(real(zVmx(nband,nband,ivk2+1,1), dp), nband = 1, ndimEV)]
   vVy(:,ivk2+1,1)=[(real(zVmy(nband,nband,ivk2+1,1), dp), nband = 1, ndimEV)]

enddo
!$omp end parallel do

ivk2=0
!$omp parallel do &
!$omp private(vk,ivk1,zh,energia,zvxy1,zvxy2,zvxyI,zevTB) &
!$omp shared(ivk2,t1,t2,t3,il,iu,coord,alpha,potential,zfock,ind_j,ind_l,numkPlot) &
!$omp shared(vEnergy,vMomentum,zVmx,zVmy)
 do ivk1 = 1,numkPlot

    vk(:)=vkBZ(ivk1+1,ivk2+1,:)

   call ComplexHk_MF(zh(:,:,nspin),coord,potential(:,nspin),alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]))

   call diagonalize(zh(:,:,nspin),energia(:,nspin),'V',il,iu)

   vEnergy(:,ivk1+1,3)=energia(:,nspin)
   vMomentum(:,ivk1+1,3)=vk(:)

   zevTB(:,il:iu)=zh(:,1:ndimEV,nspin)

   call vxy12I_MF(zvxy1,zvxy2,zvxyI,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),1_dp)

   call matrix_elements(zvxy1+zvxy2+zvxyI,zevTB,zVmx(:,:,ivk1+1,3),il,iu,il,iu)

   call vxy12I_MF(zvxy1,zvxy2,zvxyI,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),2_dp)

   call matrix_elements(zvxy1+zvxy2+zvxyI,zevTB,zVmy(:,:,ivk1+1,3),il,iu,il,iu)

enddo
!$omp end parallel do

do ivk1=1,numkPlot
write(*,*) ivk1
!$omp parallel do &
!$omp private(vk,ivk2,zh,energia,zvxy1,zvxy2,zvxyI,zevTB,nband) &
!$omp private(imu,iT,fmu,T) &
!$omp shared(ivk1,t1,t2,t3,il,iu,coord,alpha,potential,zfock,ind_j,ind_l,numkPlot,CBand,vkBZ) &
!$omp shared(vEnergy,vMomentum,zVmx,zVmy,vVDx,vVDy,vVDxNO,vVDyNO) &
!$omp reduction(+:DiamagX) &
!$omp reduction(+:DiamagY) &
!$omp reduction(+:fKinetic) &
!$omp reduction(+:DiamagNO)
 do ivk2 = 1,numkPlot

    vk(:)=vkBZ(ivk1+1,ivk2+1,:)

   call ComplexHk_MF(zh(:,:,nspin),coord,potential(:,nspin),alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]))

   call diagonalize(zh(:,:,nspin),energia(:,nspin),'V',il,iu)

   vEnergy(:,ivk2+1,2)=energia(:,nspin)
   vMomentum(:,ivk2+1,2)=vk(:)

   zevTB(:,il:iu)=zh(:,1:ndimEV,nspin)

   call vxy12I_MF(zvxy1,zvxy2,zvxyI,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),1_dp)

   call matrix_elements(zvxy1+zvxy2+zvxyI,zevTB,zVmx(:,:,ivk2+1,2),il,iu,il,iu)

   call vxy12I_MF(zvxy1,zvxy2,zvxyI,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),2_dp)

   call matrix_elements(zvxy1+zvxy2+zvxyI,zevTB,zVmy(:,:,ivk2+1,2),il,iu,il,iu)

   call vxy12D_MF(zvxy1,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),1_dp)

   call z_matrix_elements_diag_real(zvxy1,zevTB,vVDx)

   call vxy12D_MF(zvxy1,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),2_dp)

   call z_matrix_elements_diag_real(zvxy1,zevTB,vVDy)

   call vxy12ID_NO(zvxy1,zvxy2,zvxyI,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),1_dp)

   call z_matrix_elements_diag_real(zvxy1+zvxy2+zvxyI,zevTB,vVDxNO)

   call vxy12ID_NO(zvxy1,zvxy2,zvxyI,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),2_dp)

   call z_matrix_elements_diag_real(zvxy1+zvxy2+zvxyI,zevTB,vVDyNO)

   vVx(:,ivk2+1,2)=[(real(zVmx(nband,nband,ivk2+1,2), dp), nband = 1, ndimEV)]
   vVy(:,ivk2+1,2)=[(real(zVmy(nband,nband,ivk2+1,2), dp), nband = 1, ndimEV)]

   do iT=1,nT ! iT

   T=ftemperature(iT)!*0.0000861733_dp ! Temperature

   fmu=efermi

   do nband=1,ndimEV

    DiamagNO(iT)=DiamagNO(iT)+(vVDxNO(nband)+vVDyNO(nband))/2._dp*fermi_dist(energia(nband,nspin)-fmu,T)

    enddo

   do imu=1,nmu

   fmu=efermi+CBand*(nmu-imu)/real(nmu,8)

       do nband=1,ndimEV

    DiamagX(imu,iT)=DiamagX(imu,iT)+vVDx(nband)*fermi_dist(energia(nband,nspin)-fmu,T)
    DiamagY(imu,iT)=DiamagY(imu,iT)+vVDy(nband)*fermi_dist(energia(nband,nspin)-fmu,T)
    fKinetic(imu,iT)=fKinetic(imu,iT)+energia(nband,nspin)*fermi_dist(energia(nband,nspin)-fmu,T)

    enddo

    enddo
    enddo

enddo
!$omp end parallel do

vEnergy(:,1,2)=vEnergy(:,ivk1+1,3)
vMomentum(:,1,2)=vMomentum(:,ivk1+1,3)

vVx(:,1,2)=vVx(:,ivk1+1,3)

vVy(:,1,2)=vVy(:,ivk1+1,3)

zVmx(:,:,1,2)=zVmx(:,:,ivk1+1,3)

zVmy(:,:,1,2)=zVmy(:,:,ivk1+1,3)

write(*,*) 'Test'

! Set spectral range

!$omp parallel do &
!$omp private(ivk2,pivotEnergy,pivotMomentum,pivotVx,pivotVy) &
!$omp private(tempD) &
!$omp shared(numkPlot,vEnergy,vMomentum,vVx,vVy,nbounds) &
!$omp reduction(+:dos) &
!$omp reduction(+:drude)
do ivk2=0,numkPlot-1

pivotEnergy=reshape([vEnergy(:,ivk2+1,1),vEnergy(:,ivk2+2,1),vEnergy(:,ivk2+1,2),vEnergy(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVx=reshape([vVx(:,ivk2+1,1),vVx(:,ivk2+2,1),vVx(:,ivk2+1,2),vVx(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVy=reshape([vVy(:,ivk2+1,1),vVy(:,ivk2+2,1),vVy(:,ivk2+1,2),vVy(:,ivk2+2,2)],[ndimEV,4_dp])

pivotMomentum=reshape([vMomentum(:,ivk2+1,1),vMomentum(:,ivk2+2,1),vMomentum(:,ivk2+1,2),vMomentum(:,ivk2+2,2)],[2_dp,4_dp])

call dos_F(tempD,pivotMomentum,pivotEnergy,nomegaD,nspacingD,nshiftD,efermi,nbounds)

dos(:)=dos(:)+tempD(:)

call drude_F(tempD,pivotMomentum,pivotEnergy,pivotVx,pivotVy,nomegaD,nspacingD,nshiftD,efermi,nbounds)

drude(:)=drude(:)+tempD(:)

enddo ! icount
!$omp end parallel do

! Conductivity

 allocate(pivotVmx(ndimEV,ndimEV,4),pivotVmy(ndimEV,ndimEV,4))

do imu=1,nmu

fmu=efermi+CBand*(nmu-imu)/real(nmu,8)

do iT=1,nT ! iT q

T=ftemperature(iT)!*0.0000861733_dp ! Temperature

!$omp parallel do &
!$omp private(ivk2,pivotEnergy,pivotMomentum,pivotVmx,pivotVmy) &
!$omp private(temp1) &
!$omp shared(numkPlot,imu,iT,fmu,T,vEnergy,vMomentum,zVmx,zVmy,nomega,nshift,nbounds) &
!$omp reduction(+:sigmaTot)
do ivk2=0,numkPlot-1

pivotEnergy=reshape([vEnergy(:,ivk2+1,1),vEnergy(:,ivk2+2,1),vEnergy(:,ivk2+1,2),vEnergy(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVmx=reshape([zVmx(:,:,ivk2+1,1),zVmx(:,:,ivk2+2,1),zVmx(:,:,ivk2+1,2),zVmx(:,:,ivk2+2,2)],[ndimEV,ndimEV,4_dp])

pivotVmy=reshape([zVmy(:,:,ivk2+1,1),zVmy(:,:,ivk2+2,1),zVmy(:,:,ivk2+1,2),zVmy(:,:,ivk2+2,2)],[ndimEV,ndimEV,4_dp])

pivotMomentum=reshape([vMomentum(:,ivk2+1,1),vMomentum(:,ivk2+2,1),vMomentum(:,ivk2+1,2),vMomentum(:,ivk2+2,2)],[2_dp,4_dp])

call conductivityS(temp1,pivotMomentum,pivotEnergy,pivotVmx,pivotVmy,nomega,nspacing,nshift,fmu,T,nbounds)

sigmaTot(:,imu,iT)=sigmaTot(:,imu,iT)+temp1(:)

enddo ! icount
!$omp end parallel do

enddo ! iT

enddo ! imu

vEnergy(:,:,1)=vEnergy(:,:,2)
vMomentum(:,:,1)=vMomentum(:,:,2)

vVx(:,:,1)=vVx(:,:,2)

vVy(:,:,1)=vVy(:,:,2)

zVmx(:,:,:,1)=zVmx(:,:,:,2)

zVmy(:,:,:,1)=zVmy(:,:,:,2)

deallocate(pivotVmx,pivotVmy)

enddo ! ivk1

write(*,*) 'Output'

write(*,*) 'WriteDos'

write(filename,'(A,A12,I0,A9,I0,A)') trim(dir),'Dos-numkPlot',numkPlot,'-nspacing',nspacingD,trim(parameters)
open(89,file=filename,status='replace')

write(filename,'(A,A14,I0,A9,I0,A)') trim(dir),'Drude-numkPlot',numkPlot,'-nspacing',nspacingD,trim(parameters)
open(99,file=filename,status='replace')

do n=1,nomegaD
om=real(n-nshiftD)/real(nspacingD)
write(89,*) om,dos(n)/Ac,2.*abs(om)/pi/(3./4.)/2.7/2.7 ! Spin is not included
write(99,*) om,drude(n)/Ac ! Spin is not included
enddo

close(99)
close(89)

! normalization: int DOS*sqrt(3)/2*ndim/4= number of bands =iu-il

close(99)
close(79)
close(89)

deallocate(dos)

write(*,*) 'WriteSigma'

do imu=1,nmu

fmu=CBand*(nmu-imu)/real(nmu,8)

write(*,*) imu,fmu

do iT=1,nT ! iT

T=ftemperature(iT)*0.0000861733_dp

write(filename,'(A,A15,I0,A3,I0,A9,I0,A9,I0,A)') trim(dir),'ConductivityTot',int(ftemperature(iT)),'-Mu',imu,'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(87,file=filename,status='replace')

! Normalization to "sigma_0=e^2/hbar/8"
! Implementation with "exact" symmetry for negative frequencies
sumSigma=0._dp
do n=nshift,nomega
om=real(n-nshift)/real(nspacing)
write(87,*) om,sigmaTot(n,imu,iT)*8_dp/Ac
enddo
close(87)

enddo !iT

enddo !imu

write(*,*) 'Test'

write(filename,'(A,A5,A)') trim(dir),'CBand',trim(parameters)
open(87,file=filename,status='replace')

write(87,*) CBand,numkPlot,nspacing

write(87,*) DiamagNO(nT)/AreaMoire

close(87)

do iT=1,nT ! iT

T=ftemperature(iT)*0.0000861733_dp

write(filename,'(A,A6,I0,A9,I0,A9,I0,A)') trim(dir),'DiaTot',int(ftemperature(iT)),'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(67,file=filename,status='replace')

write(filename,'(A,A7,I0,A9,I0,A9,I0,A)') trim(dir),'Kinetic',int(ftemperature(iT)),'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(77,file=filename,status='replace')

do imu=1,nmu

fmu=CBand*(nmu-imu)/real(nmu,8)

write(filename,'(A,A5,I0,A3,I0,A9,I0,A9,I0,A)') trim(dir),'DiaXY',int(ftemperature(iT)),'-Mu',imu,'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(76,file=filename,status='replace')

! Normalization to "sigma_0=e^2/hbar/8"
! Implementation with "exact" symmetry for negative frequencies

sumSigma=0._dp
do n=nshift,nomega
sumSigma=sumSigma+sigmaTot(n,imu,iT)
enddo

write(67,*) fmu,-(DiamagX(imu,iT)+DiamagY(imu,iT))/2._dp/AreaMoire-sumSigma/Ac/real(nspacing,8)*(2._dp/pi)!+DiamagNO(iT)/AreaMoire

write(76,*) nomega,nshift,-DiamagX(imu,iT)/AreaMoire,-DiamagY(imu,iT)/AreaMoire,-DiamagNO(iT)/AreaMoire

write(77,*) fmu,-fKinetic(imu,iT)/AreaMoire*pi/2._dp*8_dp/6._dp,sumSigma/Ac/real(nspacing,8)*8_dp

close(76)

enddo ! iT

enddo ! imu

close(67)
close(77)

! Kramers-Kronig postprocessing: compute_kramers_kronig.f90.

write(*,*) 'TestEnd'

end program compute_spectral_optimize
