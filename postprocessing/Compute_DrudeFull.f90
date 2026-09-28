! Cleaned from plotDrudeFull.f90; original modules retained.
! Input zFock is the one-particle reduced density matrix from FORGE.
! Geometry/input: FORGE dev; response kernels: Hamiltonian.f90. See README.md.

program compute_drude_full

use omp_lib
use Hamiltonian
use PostProcessingInput, only: BuildGeometry, ReadFock, InputParameters, ValidateResponseModel
use Geometry, only: NumberNeighborCells, OrderNeighborCells
use Setup, only: alphaInput => alpha, alphaHInput => alphaH
use lapack_routines

 implicit none

 character(1024) :: filename
 character(220) :: parameters
 integer(dp) :: ivk1, ivk2, n, m, i, j

 integer(dp) :: icount, ncount, numkPlot
 integer(dp) :: ndimEV, nNP, il, iu, nband, nspin, it
 integer(dp) :: ind
 real(dp) :: fmu

 real(dp) :: vk(2), t1(2), t2(2), t3(2), vq1(2), vq12(2)
 real(dp) :: ang_mat(2,2)
 integer(dp), allocatable :: npointsBZ(:,:), npointsBack(:,:)
 integer(dp), allocatable :: ind_j(:), ind_l(:)
 real(dp), allocatable :: coord(:,:), energia(:,:), vkBZ(:,:,:)
 real(dp), allocatable :: bandsTB(:,:,:), vx1TB(:,:), vy1TB(:,:), vx2TB(:,:), vy2TB(:,:), vxITB(:,:)
 real(dp), allocatable :: vyITB(:,:)
 complex(dp), allocatable :: zh(:,:,:), zvxy1(:,:), zvxy2(:,:), zvxyI(:,:)

 real(dp), allocatable :: dos(:)
 real(dp), allocatable :: drudeTot(:), drudeMag(:), drudeChi(:)

 real(dp), allocatable :: dens(:,:)
 real(dp), allocatable :: potential(:,:)
 real(dp), allocatable :: fLongRange(:,:)
 complex(dp), allocatable :: zfock(:,:,:,:)


 real(dp) :: UPlot

 real(dp) :: alpha, alphaH

 real(dp) :: OmegaMin, OmegaMax
 integer(dp), parameter :: nspacing=100000
 integer(dp) :: nomega
 integer(dp) :: nshift

 integer(dp) :: nbounds(4)

 call ValidateResponseModel()

 ndimEV=numb
il=(ndim-2*ndimEV+2*ncb)/2+1
iu=(ndim+2*ncb)/2
nNP=(ndim/2-il)+1

 numkPlot=30

allocate(vkBZ(1:numkPlot+1,1:numkPlot+1,1:2))

allocate(dens(1:ndim,1:numS))
allocate(energia(ndim,1:numS))
allocate(coord(ndim,3))
allocate(zh(ndim,ndim,1:numS))
allocate(zvxy1(ndim,ndim),zvxy2(ndim,ndim),zvxyI(ndim,ndim))

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
allocate(bandsTB(ndimEV,ncount,1))
allocate(vx1TB(ndimEV,ncount),vx2TB(ndimEV,ncount),vxITB(ndimEV,ncount))
allocate(vy1TB(ndimEV,ncount),vy2TB(ndimEV,ncount),vyITB(ndimEV,ncount))

! Geometry: shared implementation of Main_BandStructure's construction.
call BuildGeometry(coord, t1, t2, t3, vq1, vq12, ang_mat)

! The long-range interaction uses the shared dev screening routine.
call longRange(fLongRange,coord,ndim,t1,t2)

! Get Points

nbounds(1)=1 ! min for i_v
nbounds(2)=ndimEV-ncb ! max for i_v
nbounds(3)=ndimEV-ncb+1 ! min for i_c
nbounds(4)=ndimEV ! max for i_c

! Set filename

UPlot = U
nbounds(2) = ndimEV-ncb
nbounds(3) = ndimEV-ncb+1
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

write(*,*) 'Get Points DOS'

ncount=numkPlot*numkPlot
allocate(npointsBZ(ncount,8))
allocate(npointsBack(numkPlot,numkPlot))

call samplePoints1(npointsBZ,npointsBack,numkPlot,ncount)

do ivk2=0,numkPlot-1
do ivk1=0,numkPlot-1

vkBZ(ivk1+1,ivk2+1,1)=(ivk1*vq1(1)+ivk2*vq12(1))/real(numkPlot,dp)
vkBZ(ivk1+1,ivk2+1,2)=(ivk1*vq1(2)+ivk2*vq12(2))/real(numkPlot,dp)

enddo
enddo

nspin=1
write(*,*) 'Start DOS ndimEV',ndimEV

vx1TB=0._dp
vx2TB=0._dp
vxITB=0._dp
vy1TB=0._dp
vy2TB=0._dp
vyITB=0._dp

!$omp parallel do &
!$omp private(nband,vk,icount,zh,energia,n,m,zvxy1,zvxy2,zvxyI) &
!$omp shared(ivk1,coord,numkPlot,t1,t2,t3,il,iu,bandsTB,ind_j,ind_l,ind) &
!$omp shared(potential,zfock,vx1TB,vy1TB,vx2TB,vy2TB,vxITB,vyITB,ndimEV)
 do icount = 1,ncount

    vk(:)=vkBZ(npointsBZ(icount,1)+1,npointsBZ(icount,2)+1,:)

   call ComplexHk_MF(zh(:,:,nspin),coord,potential(:,nspin),alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]))

   call diagonalize(zh(:,:,nspin),energia(:,nspin),'V',il,iu)

   do nband=1,ndimEV
   bandsTB(nband,icount,nspin)=energia(nband,nspin)
   enddo

   call vxy12I_MF(zvxy1,zvxy2,zvxyI,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),1_dp)

   call z_matrix_elements_diag_real(zvxy1,zh(:,:,nspin),vx1TB(:,icount))

   call z_matrix_elements_diag_real(zvxy2,zh(:,:,nspin),vx2TB(:,icount))

   call z_matrix_elements_diag_real(zvxyI,zh(:,:,nspin),vxITB(:,icount))

   call vxy12I_MF(zvxy1,zvxy2,zvxyI,coord,alpha,zfock(:,:,:,nspin),ind_j,ind_l,ndim,ind,vk,reshape([t1,t2,t3],[2,3]),2_dp)

   call z_matrix_elements_diag_real(zvxy1,zh(:,:,nspin),vy1TB(:,icount))

   call z_matrix_elements_diag_real(zvxy2,zh(:,:,nspin),vy2TB(:,icount))

   call z_matrix_elements_diag_real(zvxyI,zh(:,:,nspin),vyITB(:,icount))

 end do
!$omp end parallel do

write(*,*) 'End DOS'

fmu=(bandsTB(nbounds(3),npointsBack(numkPlot/3+1,numkPlot/3+1),nspin)+bandsTB(nbounds(2),npointsBack(numkPlot/3+1,numkPlot/3+1),nspin))/2.
write(*,*) fmu
write(*,*) maxval(bandsTB(ndimEV,:,nspin)),minval(bandsTB(1_dp,:,nspin))

  OmegaMax=maxval(bandsTB(ndimEV,:,nspin))-fmu
  OmegaMin=minval(bandsTB(1_dp,:,nspin))-fmu

  OmegaMax=OmegaMax*1.1_dp
  OmegaMin=OmegaMin*1.1_dp

  write(*,*) 'maxmin',OmegaMax,OmegaMin

  nomega=int((OmegaMax-OmegaMin)*nspacing)
  nshift=int(-OmegaMin*nspacing)

  write(*,*) OmegaMax,OmegaMin
  write(*,*) nomega,nshift

  allocate(dos(nomega))
  allocate(drudeTot(nomega),drudeMag(nomega),drudeChi(nomega))

write(filename,'(A,A5,A)') trim(dir),'DosI0',trim(parameters)
open(99,file=filename,status='replace')

call dosFull(bandsTB(:,:,1),dos,fmu,npointsBack,ndimEV,numkPlot,ncount,vq1,vq12,nomega,nspacing,nshift,nbounds)

write(*,*) 'Write Drude'

  write(*,*) OmegaMax,OmegaMin
  write(*,*) nomega,nshift

write(filename,'(A,A9,I0,A)') trim(dir),'DrudeTotI',numkPlot,trim(parameters)
open(99,file=filename,status='replace')

call getDrudeFull(bandsTB(:,:,1),drudeTot,drudeMag,drudeChi,vx1TB,vy1TB,vx2TB,vy2TB,vxITB,vyITB, &
fmu,npointsBack,ndimEV,numkPlot,ncount,vq1,vq12,nomega,nspacing,nshift,nbounds)

write(filename,'(A,A5,A)') trim(dir),'Bands',trim(parameters)
open(99,file=filename,status='replace')

call plotBands(bandsTB(:,:,nspin),fmu,npointsBack,ndimEV,numkPlot,ncount,vq1,vq12)

write(filename,'(A,A3,A)') trim(dir),'Vx1',trim(parameters)
open(99,file=filename,status='replace')

call plotBands(vx1TB,fmu,npointsBack,ndimEV,numkPlot,ncount,vq1,vq12)

write(filename,'(A,A3,A)') trim(dir),'Vx2',trim(parameters)
open(99,file=filename,status='replace')

call plotBands(vx2TB,fmu,npointsBack,ndimEV,numkPlot,ncount,vq1,vq12)

write(filename,'(A,A3,A)') trim(dir),'VxI',trim(parameters)
open(99,file=filename,status='replace')

call plotBands(vxITB,fmu,npointsBack,ndimEV,numkPlot,ncount,vq1,vq12)

write(filename,'(A,A9,A)') trim(dir),'BandsAllI',trim(parameters)
open(99,file=filename,status='replace')

call plotBandsAll(bandsTB(:,:,nspin),fmu,npointsBack,ndimEV,numkPlot,ncount,vq1,vq12)

close(99)

write(filename,'(A,A17,I0,A)') trim(dir),'DrudeCompare-numk',numkPlot,trim(parameters)
open(101,file=filename,status='replace')

write(filename,'(A,A12,I0,A)') trim(dir),'DosDens-numk',numkPlot,trim(parameters)
open(102,file=filename,status='replace')

     do iT=1,nT

write(filename,'(A,A9,I0,A5,I0,A)') trim(dir),'DrudeTotT',int(ftemperature(it)),'-numk',numkPlot,trim(parameters)
open(110+iT,file=filename,status='replace')

write(filename,'(A,A9,I0,A5,I0,A)') trim(dir),'DrudeMagT',int(ftemperature(it)),'-numk',numkPlot,trim(parameters)
open(120+iT,file=filename,status='replace')

write(filename,'(A,A9,I0,A5,I0,A)') trim(dir),'DrudeChiT',int(ftemperature(it)),'-numk',numkPlot,trim(parameters)
open(130+iT,file=filename,status='replace')

write(filename,'(A,A7,I0,A)') trim(dir),'DrudeGT',int(ftemperature(it)),trim(parameters)
open(140+iT,file=filename,status='replace')

write(filename,'(A,A11,I0,A)') trim(dir),'OccupationT',int(ftemperature(it)),trim(parameters)
open(150+iT,file=filename,status='replace')

enddo

call DrudeTemp(dos,drudeTot,drudeMag,drudeChi,nshift,nomega,nspacing,3_dp)


write(*,*) 'Final filling'


write(*,*) 'Final U'


      stop
      end program compute_drude_full
