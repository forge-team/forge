!
! Equivalent compute order parameters driver for the current FORGE workspace.
! This program reads an existing Fock matrix, constructs the real-space
! geometry, diagonalizes and computes the computes order-parameter quantities using the external
! OrderPlot module routines, restricted to a set of low-energy bands.
! Now, it obtains the order parameters of the occupied central 4 bands (i.e. tailored to the flat bands of TBG)
! as well as of the central numb bands (usually numb=20 for TBG).

program computeorderparams

use omp_lib
use Setup
use Geometry
use TightBinding
use HartreeFock
use OrderPlot
use lapack_routines

implicit none

character(250) :: filename
character(210) :: parameters
character(150) :: my_iomsg
character(1) :: dop

integer(dp) :: n, m, i, j, nspin, numNeighborCells, n1, n2
integer(dp) :: my_iostat, rcc
integer(dp), allocatable :: nUnitCell_1(:), nUnitCell_2(:)

real(dp) :: aMoire, cs, sn
real(dp) :: t1(2), t2(2), t3(2), tn(6,2), g1(2), g12(2), RotMatrix(2,2)
integer(dp) :: TnTonUnitCell12(0:6,2)

complex(dp) :: zinput

real(dp), allocatable :: Coords(:,:)
complex(dp), allocatable :: zFock(:,:,:,:)
integer(dp), allocatable :: NearestNeighborsUC(:,:), NearestNeighborsT(:,:)
integer(dp), allocatable :: KekuleNeighbors(:,:,:), KekuleLattice(:,:,:)
complex(dp), allocatable :: fKA_conjxfKpB(:), fKB_conjxfKpA(:), fKp_conjxfK(:)
real(dp), allocatable :: ValleyPol(:), NormSquared(:)

allocate(Coords(ndim,3))
allocate(NearestNeighborsUC(ndim,3))
allocate(NearestNeighborsT(ndim,3))
allocate(ValleyPol(ndim))
allocate(NormSquared(ndim))

call NumberNeighborCells(numI, numNeighborCells)
allocate(zFock(ndim,ndim,numNeighborCells,numS))
zFock(:,:,:,:) = cmplx(0.0_dp,0.0_dp,dp)

allocate(nUnitCell_1(numNeighborCells))
allocate(nUnitCell_2(numNeighborCells))
call OrderNeighborCells(numI, numNeighborCells, nUnitCell_1, nUnitCell_2)

! real-space geometry

aMoire = 3.0_dp*ntheta**2 + 3.0_dp*ntheta + 1.0_dp

g1  =  (4.0_dp*pi/3.0_dp)/aMoire*(real(3*ntheta+1,dp)*a1+a2)
g12 =  (4.0_dp*pi/3.0_dp)/aMoire*(real(3*ntheta+2,dp)*a2-a1)

cs = 1.0_dp-1.0_dp/(2.0_dp*aMoire)
sn = sqrt(1.0_dp-cs**2)
RotMatrix = reshape([cos(0.5_dp*acos(cs)),-sin(0.5_dp*acos(cs)),sin(0.5_dp*acos(cs)),cos(0.5_dp*acos(cs))],[2,2])

g1  = matmul(RotMatrix,g1)
g12 = matmul(RotMatrix,g12)

! Moire lattice parameters

t1 = real(ntheta,dp)*a1 + real(ntheta+1,dp)*a2
t2 = real(-ntheta-1,dp)*a1 + real(2*ntheta+1,dp)*a2
t3 = t2 - t1

call WignerSeitzCell(Coords,t1,t2,cs,sn)

do n=1,ndim
    Coords(n,:) = matmul(RotMatrix,Coords(n,1:2))
end do

t1 = matmul(RotMatrix,t1)
t2 = matmul(RotMatrix,t2)
t3 = t2 - t1

tn(1,:) = t1
tn(2,:) = t2
tn(3,:) = t3
tn(4,:) = -t1
tn(5,:) = -t2
tn(6,:) = -t3

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

call getNearestNeighbors(NearestNeighborsUC,NearestNeighborsT,Coords,ndim,tn)

! Map tn index (0-6) to (n1,n2) moiré lattice Coordsinates
TnTonUnitCell12(0,:) = [ 0,  0]
TnTonUnitCell12(1,:) = [ 1,  0]
TnTonUnitCell12(2,:) = [ 0,  1]
TnTonUnitCell12(3,:) = [-1,  1]
TnTonUnitCell12(4,:) = [-1,  0]
TnTonUnitCell12(5,:) = [ 0, -1]
TnTonUnitCell12(6,:) = [ 1, -1]

! Kekulé order-parameter geometry

allocate(KekuleNeighbors(ndim/2,2,3))
allocate(KekuleLattice(ndim/2,6,2))

call getKekuleLattice(Coords,ndim,RotMatrix,a1,a2,tn,KekuleLattice)
call getKekuleNeighbors(KekuleNeighbors,Coords,ndim,a1,a2,RotMatrix,tn)

! Read the Fock matrix from disk

if(nfilling.LT.0._dp)then
dop = '-'
else
dop = '+'
endif

write(parameters,'(A1,A3,A8,A1,F3.1,A5,I0,A2,I0,A8,I0,A6,I0,A4,F0.1,A2,F4.2,A7,I0,A3,F0.1,A6,F0.3,A5,I0,A5,I0,A5,I0,A3,I0,A4)') &
    '-',statename,'-filling',dop,abs(nfilling),'-numS',numS,'-i',ntheta,'-nlayers',nlayers,'-relax',nrelax,&
    '-eps',epsilon,'-U',U,'-screen',nscreen,'-xi',xi*0.246_dp,'-delta',Delta,&
    '-numI',numI,'-numC',numC,'-numk',numk,'-dp',dp,'.dat'

write(*,*) 'parameters', parameters

if(numS.EQ.2)then
    write(*,*) 'reading Fock for both spin components'
endif

do nspin = 1, numS
     if(numS.EQ.1)then
         write(filename,'(A9,A10,I0,A210)') dirFock,'Fock-nspin',1,parameters
     else
         write(filename,'(A9,A10,I0,A210)') dirFock,'Fock-nspin',nspin,parameters
     endif

     open(12, file=filename, form='unformatted', status='old', access='direct', recl=dp*2)
     rcc = 0_dp
     do i = 1, ndim
         do j = 1, i-1
             rcc = rcc + 1_dp
             read(12,rec=rcc) zinput
             zFock(i,j,1,nspin) = zinput
             zFock(j,i,1,nspin) = conjg(zinput)
         end do
         rcc = rcc + 1_dp
         read(12,rec=rcc) zinput
         zFock(i,i,1,nspin) = zinput
     end do
     do m = 2, numNeighborCells
         do i = 1, ndim
             do j = 1, ndim
                 rcc = rcc + 1_dp
                 read(12,rec=rcc) zinput
                 zFock(i,j,m,nspin) = zinput
             end do
         end do
     end do
     close(12)
end do

! Allocate order-parameter output arrays

allocate(fKA_conjxfKpB(ndim/2))
allocate(fKB_conjxfKpA(ndim/2))
allocate(fKp_conjxfK(ndim))

! Compute order-parameter quantities for each spin sector

do nspin = 1, numS

    call InterSubInterValAlt(fKA_conjxfKpB,fKB_conjxfKpA,KekuleLattice,KekuleNeighbors,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A21,I0,A6,I0,A210)') dir,'InterSubInterVal-numb',ndim,'-nspin',nspin,parameters
    open(98,file=filename,status='replace')
    do i = 1, ndim/2
        write(98,fmt='(E12.5,3X,E12.5,3X,E12.5,3X,E12.5)') real(fKA_conjxfKpB(i)), aimag(fKA_conjxfKpB(i)), real(fKB_conjxfKpA(i)), aimag(fKB_conjxfKpA(i))
    end do
    close(98)

    call IntraSubInterValAlt(fKp_conjxfK,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,NearestNeighborsUC,NearestNeighborsT,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A21,I0,A6,I0,A210)') dir,'IntraSubInterVal-numb',ndim,'-nspin',nspin,parameters
    open(98,file=filename,status='replace')
    do i = 1, ndim
        write(98,*) real(fKp_conjxfK(i)), aimag(fKp_conjxfK(i))
    end do
    close(98)

    call InterSubIntraValAlt(fKA_conjxfKpB,fKB_conjxfKpA,KekuleLattice,KekuleNeighbors,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A21,I0,A6,I0,A210)') dir,'InterSubIntraVal-numb',ndim,'-nspin',nspin,parameters
    open(98,file=filename,status='replace')
    do i = 1, ndim/2
        write(98,*) real(fKA_conjxfKpB(i)), aimag(fKA_conjxfKpB(i))
    end do
    do i = 1, ndim/2
        write(98,*) real(fKB_conjxfKpA(i)), aimag(fKB_conjxfKpA(i))
    end do
    close(98)

    call IntraSubIntraValAlt(ValleyPol,NormSquared,ndim,NearestNeighborsUC,NearestNeighborsT,numNeighborCells,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A21,I0,A6,I0,A210)') dir,'IntraSubIntraVal-numb',ndim,'-nspin',nspin,parameters
    open(98,file=filename,status='replace')
    do i = 1, ndim
        write(98,*) ValleyPol(i), NormSquared(i)
    end do
    close(98)

end do

end program computeorderparams
