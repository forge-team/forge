!
! Band structure along the high-symmetry path of PlotBands, on a finer k-grid than the
! self-consistent run.
! This program reads an existing Fock matrix, rebuilds the Hartree-Fock Hamiltonian and
! diagonalizes it only at the k-points of a numkPostProc x numkPostProc grid that PlotBands reads. The rest of
! Bands is never computed nor read, so numkPostProc can be much larger than the numk of the Fock matrix.
! The Fermi energy is the last value in the Mu output of the self-consistent run.

program compute_plotBands

use omp_lib
use Setup
use Geometry
use TightBinding
use HartreeFock
use lapack_routines
use PostProcessingInput, only: BuildGeometry, InputParameters, ReadFock

implicit none


! Plot grid and path options
integer(dp), parameter :: numkPostProc = 60
character(240) :: outputSuffix
integer(dp), parameter :: nMirrorPath = 0    ! 0: Kp -> G -> K -> M -> G -> K / 1: mirrored path  K -> G -> Kp -> M -> G -> Kp

character(300) :: filename
character(210) :: inputSuffix
integer(dp), parameter :: NkB = numkPostProc*numkPostProc
integer(dp) :: n, m, i, j, nspin, numNeighborCells, rcc, my_iostat
integer(dp) :: ivk, ivk1, ivk2, jvk1, jvk2, iflat, ipath, ipoint, numPath, numPoints

real(dp) :: FermiEnergy, itinput, muinput
real(dp) :: t1(2), t2(2), t3(2), tn(6,2), tnHF(2,3), g1(2), g12(2), RotMatrix(2,2), vk(2)


integer(dp), allocatable :: nUnitCell_1(:), nUnitCell_2(:)
integer(dp), allocatable :: NearestNeighborsUC(:,:), NearestNeighborsT(:,:)
integer(dp), allocatable :: nMomentaComponents(:,:), nMomentaFlattened(:,:)
integer(dp), allocatable :: nPath(:,:), nPoints(:,:)
logical, allocatable :: lOnPath(:)

real(dp), allocatable :: Coords(:,:), MomentaValues(:,:,:)
real(dp), allocatable :: Density(:,:), DensitySub(:), LongRange(:,:), Potential(:,:)
real(dp), allocatable :: Bands(:,:,:), Energies(:)

complex(dp), allocatable :: zFock(:,:,:,:), zH(:,:)

if(mod(numkPostProc,6_dp).NE.0)then
    write(*,*) 'ERROR: numkPostProc must be a multiple of 6, numkPostProc = ', numkPostProc
    stop 1
endif

allocate(Coords(ndim,3))
allocate(NearestNeighborsUC(ndim,3))
allocate(NearestNeighborsT(ndim,3))
allocate(Density(ndim,numS))
allocate(DensitySub(ndim))
allocate(Potential(ndim,numS))
allocate(LongRange(ndim,ndim))

allocate(nMomentaComponents(NkB,2))
allocate(nMomentaFlattened(numkPostProc,numkPostProc))
allocate(MomentaValues(numkPostProc,numkPostProc,2))
allocate(Bands(numb,NkB,numS))
Bands(:,:,:) = 0.0_dp

call NumberNeighborCells(numI, numNeighborCells)
allocate(zFock(ndim,ndim,numNeighborCells,numS))
zFock(:,:,:,:) = cmplx(0.0_dp,0.0_dp,dp)

allocate(nUnitCell_1(numNeighborCells))
allocate(nUnitCell_2(numNeighborCells))
call OrderNeighborCells(numI, numNeighborCells, nUnitCell_1, nUnitCell_2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Geometry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! reciprocal lattice vectors of the superlattice
call BuildGeometry(Coords,t1,t2,t3,g1,g12,RotMatrix)
call SampleBZ(nMomentaComponents,nMomentaFlattened,MomentaValues,numkPostProc,g1,g12)

tn(1,:) = t1
tn(2,:) = t2
tn(3,:) = t3
tn(4,:) = -t1
tn(5,:) = -t2
tn(6,:) = -t3
call getNearestNeighbors(NearestNeighborsUC,NearestNeighborsT,Coords,ndim,tn)

tnHF = reshape([t1,t2,t3],[2,3])

call LongRangeInteraction(LongRange,Coords,ndim,t1,t2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read the Fock matrix from disk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call InputParameters(inputSuffix)
write(*,*) 'parameters', trim(inputSuffix)
call ReadFock(zFock,numNeighborCells,inputSuffix)
write(outputSuffix,'(A,I0,A)') '-numkPostProc',numkPostProc,trim(inputSuffix)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Hartree and Hubbard potential of the Fock matrix, as in Main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DensitySub(:) = 0.5_dp

do n=1,ndim
    do nspin=1,numS
        Density(n,nspin) = real(zFock(n,n,1,nspin),dp)
    enddo
enddo

Potential(:,:) = 0.0_dp
if(numS.EQ.1)then
    do n=1,ndim
    do m=1,ndim
        Potential(n,1) = Potential(n,1) + 2.0_dp*(Density(m,1)-DensitySub(m))*alphaH*LongRange(n,m)
    enddo
    enddo
    Potential(:,1) = Potential(:,1) + U*(Density(:,1)-DensitySub(:))
endif

if(numS.EQ.2)then
    do n=1,ndim
        do m=1,ndim
            Potential(n,1) = Potential(n,1) + (Density(m,1)+Density(m,2)-2.0_dp*DensitySub(m))*alphaH*LongRange(n,m)
        enddo
    enddo

    do n=1,ndim
        Potential(n,2) = Potential(n,1) + U*(Density(n,1)-DensitySub(n))
        Potential(n,1) = Potential(n,1) + U*(Density(n,2)-DensitySub(n))
    enddo
endif

deallocate(LongRange)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Fermi energy: last value in the Mu output of the self-consistent run
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FermiEnergy = 0.0_dp
write(filename,'(A,A,A)') trim(dir),'Mu',trim(inputSuffix)
open(13, file=filename, status='old', iostat=my_iostat)
if(my_iostat.EQ.0)then
    do
        read(13,*,iostat=my_iostat) itinput, muinput
        if(my_iostat.NE.0) exit
        FermiEnergy = muinput
    enddo
    close(13)
else
    write(*,*) 'WARNING: ', trim(filename), ' not found, bands are not shifted by the Fermi energy'
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Grid indices (ivk1,ivk2) read by PlotBands, in its order. Must match PlotBands (p=0).
!! With k = (ivk1*g1 + ivk2*g12)/numkPostProc: K = (g1+g12)/3, Kp = -K, M = (g1+g12)/2.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

numPath = 5*numkPostProc/3 + 1
allocate(nPath(numPath,2))

ipath = 0
do ivk=0,numkPostProc/3-1                                      ! Kp -> G
    ipath = ipath + 1
    nPath(ipath,:) = [modulo(2*numkPostProc/3+ivk,numkPostProc), modulo(2*numkPostProc/3-2*ivk,numkPostProc)]
enddo
ipath = ipath + 1                                       ! G
nPath(ipath,:) = [0_dp, 0_dp]
do ivk=1,numkPostProc/3                                        ! G -> K
    ipath = ipath + 1
    nPath(ipath,:) = [modulo(ivk,numkPostProc), modulo(numkPostProc-2*ivk,numkPostProc)]
enddo
do ivk=1,numkPostProc/6                                        ! K -> M
    ipath = ipath + 1
    nPath(ipath,:) = [modulo(numkPostProc/3+ivk,numkPostProc), modulo(numkPostProc/3+ivk,numkPostProc)]
enddo
do ivk=1,numkPostProc/2-1                                      ! M -> G
    ipath = ipath + 1
    nPath(ipath,:) = [modulo(numkPostProc/2-ivk,numkPostProc), modulo(numkPostProc/2+ivk,numkPostProc)]
enddo
do ivk=0,numkPostProc/3                                        ! G -> K
    ipath = ipath + 1
    nPath(ipath,:) = [modulo(ivk,numkPostProc), modulo(ivk,numkPostProc)]
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Distinct grid points of the path (G and K appear twice), each diagonalized once
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(lOnPath(NkB))
allocate(nPoints(numPath,2))
lOnPath(:) = .false.
numPoints = 0
do ipath=1,numPath
    iflat = nMomentaFlattened(nPath(ipath,1)+1,nPath(ipath,2)+1)
    if(.not.lOnPath(iflat))then
        lOnPath(iflat) = .true.
        numPoints = numPoints + 1
        nPoints(numPoints,:) = nPath(ipath,:)
    endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Diagonalize at the path points. Bands is filled at the index PlotBands reads; for
!! nMirrorPath=1 the Hamiltonian is evaluated at the mirrored momentum (-ivk2,-ivk1).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do nspin=1,numS

    !$omp parallel do &
    !$omp private(ipoint,ivk1,ivk2,jvk1,jvk2,vk,zH,Energies) &
    !$omp shared(nspin,numPoints,nPoints,nMomentaFlattened,MomentaValues,Coords,Potential,tnHF) &
    !$omp shared(zFock,nUnitCell_1,nUnitCell_2,numNeighborCells,NearestNeighborsUC,NearestNeighborsT,Bands)
    do ipoint=1,numPoints

        allocate(zH(ndim,ndim))
        allocate(Energies(ndim))

        ivk1 = nPoints(ipoint,1)
        ivk2 = nPoints(ipoint,2)

        if(nMirrorPath.EQ.1)then
            jvk1 = modulo(-ivk2,numkPostProc)
            jvk2 = modulo(-ivk1,numkPostProc)
        else
            jvk1 = ivk1
            jvk2 = ivk2
        endif
        vk(:) = MomentaValues(jvk1+1,jvk2+1,:)

        call HamiltonianHartreeFock(zH,Coords,Potential(:,nspin),alpha,Delta,zFock(:,:,:,nspin),nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,vk,tnHF,NearestNeighborsUC,NearestNeighborsT)

        call diagonalize(zH,Energies,'N',nLower,nUpper)

        Bands(1:numb,nMomentaFlattened(ivk1+1,ivk2+1),nspin) = Energies(1:numb)

        deallocate(zH)
        deallocate(Energies)

    enddo
    !$omp end parallel do

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Output bands
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do nspin=1,numS
    write(filename,'(A,A,I0,A)') trim(dir),'Bands-nspin',nspin,trim(outputSuffix)
    open(99,file=filename,status='replace')
    call PlotBands(Bands(:,:,nspin),FermiEnergy,nMomentaFlattened,numb,numkPostProc,NkB,g1,g12)
    close(99)
enddo

end program compute_plotBands
