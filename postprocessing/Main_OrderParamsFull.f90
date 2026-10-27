! Order parameters of the complete saved Fock matrix, using OrderPlot.

program computeorderparams

use omp_lib
use Setup
use Geometry
use TightBinding
use HartreeFock
use OrderPlot
use lapack_routines
use PostProcessingInput, only: BuildGeometry, InputParameters, ReadFock

implicit none

character(300) :: filename
character(210) :: inputSuffix
character(150) :: my_iomsg

integer(dp) :: n, m, i, j, nspin, numNeighborCells, n1, n2
integer(dp) :: my_iostat, rcc
integer(dp), allocatable :: nUnitCell_1(:), nUnitCell_2(:)

real(dp) :: t1(2), t2(2), t3(2), tn(6,2), g1(2), g12(2), RotMatrix(2,2)
integer(dp) :: TnTonUnitCell12(0:6,2)


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

call BuildGeometry(Coords,t1,t2,t3,g1,g12,RotMatrix)

tn(1,:) = t1; tn(2,:) = t2; tn(3,:) = t3
tn(4,:) = -t1; tn(5,:) = -t2; tn(6,:) = -t3

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

call InputParameters(inputSuffix)
write(*,*) 'parameters', trim(inputSuffix)
call ReadFock(zFock,numNeighborCells,inputSuffix)

! Allocate order-parameter output arrays

allocate(fKA_conjxfKpB(ndim/2))
allocate(fKB_conjxfKpA(ndim/2))
allocate(fKp_conjxfK(ndim))

! Compute order-parameter quantities for each spin sector

do nspin = 1, numS

    call InterSubInterValAlt(fKA_conjxfKpB,fKB_conjxfKpA,KekuleLattice,KekuleNeighbors,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A21,I0,A6,I0,A210)') dir,'InterSubInterVal-numb',ndim,'-nspin',nspin,inputSuffix
    open(98,file=filename,status='replace')
    do i = 1, ndim/2
        write(98,'(4(ES12.5,3X))') real(fKA_conjxfKpB(i)), aimag(fKA_conjxfKpB(i)), real(fKB_conjxfKpA(i)), aimag(fKB_conjxfKpA(i))
    end do
    close(98)

    call IntraSubInterValAlt(fKp_conjxfK,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,NearestNeighborsUC,NearestNeighborsT,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A21,I0,A6,I0,A210)') dir,'IntraSubInterVal-numb',ndim,'-nspin',nspin,inputSuffix
    open(98,file=filename,status='replace')
    do i = 1, ndim
        write(98,'(2(ES12.5,3X))') real(fKp_conjxfK(i)), aimag(fKp_conjxfK(i))
    end do
    close(98)

    call InterSubIntraValAlt(fKA_conjxfKpB,fKB_conjxfKpA,KekuleLattice,KekuleNeighbors,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A21,I0,A6,I0,A210)') dir,'InterSubIntraVal-numb',ndim,'-nspin',nspin,inputSuffix
    open(98,file=filename,status='replace')
    do i = 1, ndim/2
        write(98,'(2(ES12.5,3X))') real(fKA_conjxfKpB(i)), aimag(fKA_conjxfKpB(i))
    end do
    do i = 1, ndim/2
        write(98,'(2(ES12.5,3X))') real(fKB_conjxfKpA(i)), aimag(fKB_conjxfKpA(i))
    end do
    close(98)

    call IntraSubIntraValAlt(ValleyPol,NormSquared,ndim,NearestNeighborsUC,NearestNeighborsT,numNeighborCells,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A21,I0,A6,I0,A210)') dir,'IntraSubIntraVal-numb',ndim,'-nspin',nspin,inputSuffix
    open(98,file=filename,status='replace')
    do i = 1, ndim
        write(98,'(2(ES12.5,3X))') ValleyPol(i), NormSquared(i)
    end do
    close(98)

end do

end program computeorderparams
