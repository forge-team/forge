! ============================================================
! Program: FORGE
! Authors: Tobias Stauber (lead developer)
!          Miguel Sánchez Sánchez (co-lead developer)
!          Igor Vasilevskiy (contributor)
!          José Gonzázlez (contributor)
!          José Carlos Mouriño Gallego (contributor)
!          Martin Wackerl (contributor)
!          Paul Wenk (contributor)
!	   John Schliemann (contributor)
!           
! License: MIT
! Description:
!   FORGE is a scalable, self-consistent electronic structure toolkit designed for large-scale tight-binding models.
!   It enables Hartree-Fock-level calculations in systems with very large or non-trivial unit cells using
!   real-space methods, making it ideal for nanoscale quantum materials.
!
! Usage:
!   Compile with: ifort -qopenmp -o forge.out LapackRoutines.f90 Setup.f90 Geometry.f90 TightBinding.f90
!                 HartreeFock.f90 Main.f90 -qmkl -lpthread -lm
!
!   Set up: ulimit -s unlimited
!           export OMP_NUM_THREADS=N
!           export OMP_STACKSIZE=Xg ( for total memory X*N GByte )
!
!   Run with:     ./forge.out
!
! Repository:
!   https://github.com/forge-team/
!
! Last updated: 2025-10-10
! ============================================================

program FORGE

use omp_lib
use Setup
use Geometry
use TightBinding
use HartreeFock
use lapack_routines

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  Variable declaration !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(250) :: filename
character(210) :: parameters, parametersIn
character(150) :: my_iomsg
character(1) :: dop

integer(dp) :: Nk = numk*numk
integer(dp) :: icount, n, m, i, j, nspin, it, numNeighborCells
integer(dp) :: my_iostat, rcc
integer(dp) :: nEnergyFile,nConvergenceFile
integer(dp) :: nFermiLevel, nOccStates(numS), nPartOccStates(numS), nWindowSort_2Spins(2,2)

real(dp) :: aMoire, cs, sn
real(dp) :: t1(2), t2(2), t3(2), g1(2), g12(2), RotMatrix(2,2)
real(dp) :: DensityConvergence(2), DegFactor, FermiEnergy, Step
real(dp) :: KineticEnergy(numS),FockEnergy(numS),HartreeEnergy,HubbardEnergy,TotalEnergy
real(dp) :: KineticEnergyIn(numS),FockEnergyIn(numS),HartreeEnergyIn,HubbardEnergyIn,TotalEnergyIn

complex(dpIn) :: zinput

integer(dp), allocatable :: nMomentaComponents(:,:), nMomentaFlattened(:,:), nSortedMomenta(:,:,:)
integer(dp), allocatable :: nUnitCell_1(:), nUnitCell_2(:), nC2pairs(:), nC3pairs(:)

real(dp), allocatable :: Coords(:,:), MomentaValues(:,:,:)
real(dp), allocatable :: Bands(:,:,:), SortedEnergies_1Spin(:,:), SortedEnergies_2Spins(:,:)
real(dp), allocatable :: Density(:,:), DensityIn(:,:),DensitySub(:)
real(dp), allocatable :: LongRange(:,:), Potential(:,:)

complex(dp), allocatable :: zFockBulk(:,:,:,:), zFock(:,:,:,:), zFockIn(:,:,:,:)
complex(dp), allocatable :: zEigenvectors(:,:,:,:), zSortedEigenvectors(:,:,:)

allocate(nMomentaComponents(numk*numk,2))
allocate(nMomentaFlattened(numk,numk))
allocate(MomentaValues(1:numk,1:numk,1:2))
allocate(Coords(ndim,3))
allocate(nC2pairs(ndim))
allocate(nC3pairs(ndim))
allocate(Density(1:ndim,1:numS))
allocate(DensityIn(ndim,1:numS))
allocate(DensitySub(ndim))
allocate(Potential(1:ndim,1:numS))
Potential(:,:) = 0.0_dp
allocate(LongRange(1:ndim,1:ndim))

allocate(zEigenvectors(1:ndim,1:numb,1:numk*numk,1:numS))
allocate(Bands(1:numb,1:numk*numk,1:numS))
allocate(zSortedEigenvectors(1:ndim,1:numb*numk*numk,1:numS))
allocate(nSortedMomenta(1:numb*numk*numk,1:2,1:numS))
allocate(SortedEnergies_1Spin(1:numb*numk*numk,1:numS))

if(numS.EQ.2)then
    allocate(SortedEnergies_2Spins(1:2*numb*numk*numk,1:2))
endif

call NumberNeighborCells(numI,numNeighborCells) !set numNeighborCells

allocate(zFockBulk(1:ndim,1:ndim,1:numNeighborCells,1:numS))
allocate(zFock(1:ndim,1:ndim,1:numNeighborCells,1:numS))
allocate(zFockIn(1:ndim,1:ndim,1:numNeighborCells,1:numS))
zFockBulk(:,:,:,:) = cmplx(0.0_dp,0.0_dp,dp)
zFock(:,:,:,:) = cmplx(0.0_dp,0.0_dp,dp)
zFockIn(:,:,:,:) = cmplx(0.0_dp,0.0_dp,dp)

Density(:,:) = 0.0_dp
DensityIn(:,:) = 0.0_dp

DensitySub(:) = 0.5_dp

allocate(nUnitCell_1(1:numNeighborCells))
allocate(nUnitCell_2(1:numNeighborCells))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Geometry !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call OrderNeighborCells(numI, numNeighborCells, nUnitCell_1, nUnitCell_2)

! reciprocal lattice vectors of the superlattice
aMoire = 3.0_dp*ntheta**2 + 3.0_dp*ntheta + 1.0_dp
g1  =  (4.0_dp*pi/3.0_dp)/aMoire*(real(3*ntheta+1,dp)*a1+a2)
g12  =  (4.0_dp*pi/3.0_dp)/aMoire*(real(3*ntheta+2,dp)*a2-a1)    ! g12 =  g1 + g2 (g_i·t_j=2pi delta_ij see t_1,2 below)

! angle of rotation cs
cs = 1.0_dp-1.0_dp/(2.0_dp*aMoire)
sn = sqrt(1.0_dp-cs**2)
RotMatrix = reshape([cos(0.5_dp*acos(cs)),-sin(0.5_dp*acos(cs)),sin(0.5_dp*acos(cs)),cos(0.5_dp*acos(cs))],[2,2])

g1  = matmul(RotMatrix,g1)
g12  = matmul(RotMatrix,g12)

call SampleBZ(nMomentaComponents,nMomentaFlattened,MomentaValues,numk,g1,g12)

! Moire lattice parameters

t1 = real( ntheta  ,dp)*a1 + real(  ntheta+1,dp)*a2
t2 = real(-ntheta-1,dp)*a1 + real(2*ntheta+1,dp)*a2
t3 = t2-t1

call WignerSeitzCell(Coords,t1,t2,cs,sn)
! Rotate cell to symmetrize around x-axis
do n=1,ndim
    Coords(n,:) = matmul(RotMatrix,Coords(n,1:2))
end do

t1 = matmul(RotMatrix,t1)
t2 = matmul(RotMatrix,t2)
t3 = t2-t1

if(nrelax.EQ.1)then
    call LatticeRelaxation(Coords,g1,g12)
endif

call C2_RelatedPoints(nC2pairs,Coords,ndim)
call C3_RelatedPoints(nC3pairs,Coords,ndim)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Define interaction and the neighbourhood 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call LongRangeInteraction(LongRange,Coords,ndim,t1,t2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  Initialization !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(nfilling.LT.0._dp)then
    dop='-'
else
    dop='+'
endif

write(parameters,'(A1,A3,A8,A1,F3.1,A5,I0,A2,I0,A8,I0,A6,I0,A4,F0.1,A2,F4.2,A7,I0,A3,F0.1,A6,F0.3,A5,I0,A5,I0,A5,I0,A3,I0,A4)')&
     '-',statename,'-filling',dop,abs(nfilling),'-numS',numS,'-i',ntheta,'-nlayers',nlayers,'-relax',nrelax,&
     '-eps',epsilon,'-U',U,'-screen',nscreen,'-xi',xi*0.246_dp,'-delta',Delta,&
     '-numI',numI,'-numC',numC,'-numk',numk,'-dp',dp,'.dat'
write(*,*) 'parameters',parameters


if(nRead.EQ.1)then

    if(nfillingIn.LT.0._dp)then
        dop='-'
    else
        dop='+'
    endif
    
    write(parametersIn,'(A1,A3,A8,A1,F3.1,A5,I0,A2,I0,A8,I0,A6,I0,A4,F0.1,A2,F4.2,A7,I0,A3,F0.1,A6,F0.3,A5,I0,A5,I0,A5,I0,A3,I0,A4)')&
     '-',statenameIn,'-filling',dop,abs(nfillingIn),'-numS',numSIn,'-i',ntheta,'-nlayers',nlayers,'-relax',nrelaxIn,&
     '-eps',epsilonIn,'-U',UIn,'-screen',nscreenIn,'-xi',xiIn*0.246_dp,'-delta',DeltaIn,&
     '-numI',numIIn,'-numC',numCIn,'-numk',numkIn,'-dp',dpIn,'.dat'
    
    write(*,*) 'parametersIn',parametersIn

    do nspin=1,numS
        write(*,*) 'reading Fock...'
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(numSIn.eq.1_dp)then
            write(filename,'(A9,A10,I0,A210)') dirFock,'Fock-nspin',1,parametersIn
        else
            write(filename,'(A9,A10,I0,A210)') dirFock,'Fock-nspin',nspin,parametersIn
        endif
        open(12, file=filename , form='unformatted' , status='old',access='direct',recl=dpIn*2)
        rcc=0
        do i=1,ndim
            do j=1,i-1
            !do j=1,ndim
                rcc = rcc+1
                read(12,rec=rcc) zinput   
                zFock(i,j,1,nspin) = zinput
                zFock(j,i,1,nspin) = conjg(zinput)
            enddo
            rcc = rcc+1
            read(12,rec=rcc) zinput   
            zFock(i,i,1,nspin) = zinput
        enddo
        do m=2,numNeighborCells
            do i=1,ndim
                do j=1,ndim
                    rcc = rcc+1
                    read(12,rec=rcc) zinput   
                    zFock(i,j,m,nspin) = zinput
                enddo
            enddo
        enddo 
        close(12)        
    enddo

else

    write(*,*) 'first diagonalization...'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Enforce three-fold symmetry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(nenforceC3.eq.1)then
        call Solve_C3(zFockBulk(:,:,:,1),zFock(:,:,:,1),Potential(:,1),alpha,Bands(:,:,1), &
            zEigenvectors(:,:,:,1),Coords,MomentaValues,nMomentaComponents,nLower,nUpper,ndim,RotateLayers,numk,Nk,numNeighborCells,nUnitCell_1,nUnitCell_2,reshape([t1,t2,t3],[2,3]),g1,g12,nC3pairs)
    else
        call Solve(zFockBulk(:,:,:,1),zFock(:,:,:,1),Potential(:,1),alpha,Bands(:,:,1), &
            zEigenvectors(:,:,:,1),Coords,MomentaValues,nMomentaComponents,nLower,nUpper,ndim,numk,Nk,numNeighborCells,nUnitCell_1,nUnitCell_2,reshape([t1,t2,t3],[2,3]))
    endif

    write(*,*) 'sorting...'
    nspin=1
    call SortEnergies_1Spin(ndim,numk,numb,Nk,nMomentaComponents,zSortedEigenvectors(:,:,nspin),SortedEnergies_1Spin(:,nspin),nSortedMomenta(:,:,nspin),zEigenvectors(:,:,:,nspin),Bands(:,:,nspin))
    nFermiLevel = nint(real(numk*numk*NeutralityPoint,dp) + real(numk*numk,dp)/2.0_dp*real(nfilling,dp))
    FermiEnergy=(SortedEnergies_1Spin(nFermiLevel,1)+SortedEnergies_1Spin(nFermiLevel+1,1))/2.
    nOccStates(1)=0
    nPartOccStates(1)=0
    do icount=1,Nk*numb
        if((SortedEnergies_1Spin(icount,1) - FermiEnergy).lt.-EnergyTolerance)then
            nOccStates(1) = nOccStates(1) + 1
        else if(abs(SortedEnergies_1Spin(icount,1) - FermiEnergy).lt.EnergyTolerance)then
            nPartOccStates(1) = nPartOccStates(1) + 1
            nOccStates(1) = nOccStates(1) + 1
        endif
    enddo
    DegFactor = 1.0_dp/real(nPartOccStates(1),dp)*real(nFermiLevel+nPartOccStates(1)-nOccStates(1),dp)

    write(*,*) 'Sorted energies',nspin,nOccStates(1),nPartOccStates(1),DegFactor

    write(*,*) 'get zFock...'
    zFock(:,:,:,nspin) = cmplx(0.0_dp,0.0_dp,dp)
    if(nenforceC3.eq.1)then
        call GetFock_C3(zFock(:,:,:,nspin),zSortedEigenvectors(:,:,nspin),ndim,RotateLayers,numb,numk,&
        DegFactor,nOccStates(nspin),nPartOccStates(nspin),Coords,MomentaValues,numNeighborCells,nUnitCell_1,nUnitCell_2,nSortedMomenta(:,:,nspin),nC3pairs)
    else
        call GetFock(zFock(:,:,:,nspin),zSortedEigenvectors(:,:,nspin),ndim,numb,numk,&
        DegFactor,nOccStates(nspin),nPartOccStates(nspin),Coords,MomentaValues,numNeighborCells,nUnitCell_1,nUnitCell_2,nSortedMomenta(:,:,nspin))
    endif

    zFock(:,:,:,nspin)=zFock(:,:,:,nspin)+zFockBulk(:,:,:,nspin)
    zFock(:,:,:,nspin)=zFock(:,:,:,nspin)/cmplx(real(numk*numk,dp),0.0_dp,dp)
   
    if(numS.eq.2)then
        zFock(:,:,:,2) =  zFock(:,:,:,1)
    endif

endif

do n=1,ndim
    do nspin=1,numS
        Density(n,nspin)=real(zFock(n,n,1,nspin),dp)
    enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!   Energy  !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) 'get energies, ...'

do nspin=1,numS
        call GetKineticEnergy(KineticEnergy(nspin),Coords,zFock(:,:,:,nspin),Delta,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,reshape([t1,t2,t3],[2,3]))
        call GetFockEnergy(FockEnergy(nspin),Coords,zFock(:,:,:,nspin),nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,reshape([t1,t2,t3],[2,3]))
enddo

if(numS.EQ.1)then
    call GetHartreeHubbardEnergy(HartreeEnergy,HubbardEnergy,LongRange,Density(:,1)-DensitySub(:),Density(:,1)-DensitySub(:),ndim)
else
    call GetHartreeHubbardEnergy(HartreeEnergy,HubbardEnergy,LongRange,Density(:,1)-DensitySub(:),Density(:,2)-DensitySub(:),ndim)
endif

TotalEnergy = alpha*sum(FockEnergy)*real(3-numS,dp) + sum(KineticEnergy)*real(3-numS,dp) + U*HubbardEnergy + alphaH*HartreeEnergy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!  Potential !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Potential(:,:)=0.0_dp
if(numS.EQ.1)then
    do n=1,ndim
    do m=1,ndim
        Potential(n,1)=Potential(n,1)+2.0_dp*(Density(m,1)-DensitySub(m))*alphaH*LongRange(n,m)
    enddo
    enddo
    Potential(:,1)=Potential(:,1)+U*(Density(:,1)-DensitySub(:))
endif

if(numS.EQ.2)then

    do n=1,ndim
        do m=1,ndim
            Potential(n,1)=Potential(n,1)+(Density(m,1)+Density(m,2)-2.0_dp*DensitySub(m))*alphaH*LongRange(n,m)
        enddo
    enddo

    do n=1,ndim
        Potential(n,2)=Potential(n,1)+U*(Density(n,1)-DensitySub(n))
        Potential(n,1)=Potential(n,1)+U*(Density(n,2)-DensitySub(n))
    enddo

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DensityIn = Density
zFockIn = zFock

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! output files !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(filename,'(A7,A6,A210)') dir,'Energy',parameters
nEnergyFile=66
open(nEnergyFile+1,file=filename,status='replace')

write(filename,'(A7,A11,A210)') dir,'Convergence',parameters
nConvergenceFile=68
open(nConvergenceFile+1,file=filename,status='replace')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! Self-consistency loop !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

it = 0
do while(it.LT.itmax)
    it=it+1
    write(*,*) 'it=',it

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!! Get new fock0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) 'get zFockBulk...'
    do nspin=1,numS
    
        if(nenforceC3.eq.1)then
            call Solve_C3(zFockBulk(:,:,:,nspin),zFock(:,:,:,nspin),Potential(:,nspin),alpha,Bands(:,:,nspin), &
                zEigenvectors(:,:,:,nspin),Coords,MomentaValues,nMomentaComponents,nLower,nUpper,ndim,RotateLayers,numk,Nk,numNeighborCells,nUnitCell_1,nUnitCell_2,reshape([t1,t2,t3],[2,3]),g1,g12,nC3pairs)
        else
            call Solve(zFockBulk(:,:,:,nspin),zFock(:,:,:,nspin),Potential(:,nspin),alpha,Bands(:,:,nspin), &
                zEigenvectors(:,:,:,nspin),Coords,MomentaValues,nMomentaComponents,nLower,nUpper,ndim,numk,Nk,numNeighborCells,nUnitCell_1,nUnitCell_2,reshape([t1,t2,t3],[2,3]))
        endif


    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!! Sort energies !!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) 'Sorting...'

    do nspin=1,numS
        call SortEnergies_1Spin(ndim,numk,numb,Nk,nMomentaComponents,zSortedEigenvectors(:,:,nspin),&
        SortedEnergies_1Spin(:,nspin),nSortedMomenta(:,:,nspin),zEigenvectors(:,:,:,nspin),Bands(:,:,nspin))
    enddo

    if(numS.EQ.1)then

        nspin=1

        nFermiLevel=numS*nint(real(numk*numk*NeutralityPoint,dp)+real(numk*numk,dp)/2.0_dp*nfilling)
        FermiEnergy=(SortedEnergies_1Spin(nFermiLevel,1)+SortedEnergies_1Spin(nFermiLevel+1,1))/2.
        nOccStates(1) = 0
        nPartOccStates(1) = 0
        do icount=1,Nk*numb
            if((SortedEnergies_1Spin(icount,1) - FermiEnergy).lt.-EnergyTolerance)then
                nOccStates(1) = nOccStates(1) + 1
            else if(abs(SortedEnergies_1Spin(icount,1) - FermiEnergy).lt.EnergyTolerance)then
                nOccStates(1) = nOccStates(1) + 1
                nPartOccStates(1) = nPartOccStates(1) + 1
            endif
        enddo
        DegFactor = 1.0_dp/real(nPartOccStates(1),dp)*real(nFermiLevel+nPartOccStates(1)-nOccStates(1),dp)

        write(*,*) 'Sorted energies',nspin,nOccStates(1),nPartOccStates(1),DegFactor

    else

        nWindowSort_2Spins(1,1) = NeutralityPoint - 1
        nWindowSort_2Spins(2,1) = NeutralityPoint + 2 +1
        nWindowSort_2Spins(1,2) = NeutralityPoint - 1
        nWindowSort_2Spins(2,2) = NeutralityPoint + 2 +1

        call SortEnergies_2Spins(numk,numb,nWindowSort_2Spins,SortedEnergies_1Spin,SortedEnergies_2Spins)
        
        nFermiLevel = numS*nint(real(numk*numk*NeutralityPoint,dp)+real(numk*numk,dp)/2.0_dp*nfilling)
        FermiEnergy = (SortedEnergies_2Spins(nFermiLevel,1)+SortedEnergies_2Spins(nFermiLevel+1,1))/2.

        nOccStates(1) = 0
        nPartOccStates(1) = 0
        nOccStates(2) = 0
        nPartOccStates(2) = 0
        do icount = 1, 2*numb*Nk
            if(SortedEnergies_2Spins(icount,2).EQ.1)then
                if((SortedEnergies_2Spins(icount,1) - FermiEnergy).lt.-EnergyTolerance)then
                    nOccStates(1)  = nOccStates(1)  + 1
                else if(abs(SortedEnergies_2Spins(icount,1) - FermiEnergy).lt.EnergyTolerance)then
                    nPartOccStates(1) = nPartOccStates(1) + 1
                    nOccStates(1)  = nOccStates(1)  + 1
                endif
            else
                if((SortedEnergies_2Spins(icount,1) - FermiEnergy).lt.-EnergyTolerance)then
                    nOccStates(2)  = nOccStates(2)  + 1
                else if(abs(SortedEnergies_2Spins(icount,1) - FermiEnergy).lt.EnergyTolerance)then
                    nPartOccStates(2) = nPartOccStates(2) + 1
                    nOccStates(2)  = nOccStates(2)  + 1
                endif
            endif
        enddo

        if((nPartOccStates(1).ne.0).or.(nPartOccStates(2).ne.0))then
            DegFactor = 1.0_dp/real(nPartOccStates(1)+nPartOccStates(2),dp)*real(nFermiLevel+nPartOccStates(1)+nPartOccStates(2)-nOccStates(1)-nOccStates(2),dp)
        else
            DegFactor = 0.0_dp
        endif

        write(*,*) 'Sorted energies', 1, nOccStates(1), nPartOccStates(1), 2, nOccStates(2), nPartOccStates(2), DegFactor

    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!! Get Fock !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    do nspin=1,numS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Force Valley Symmetry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(nenforceValley.eq.1)then

            write(*,*) 'get Fock...'
            zFock(:,:,:,nspin) = cmplx(0.0_dp,0.0_dp,dp)
            if(nenforceC3.eq.1)then
                call GetFock_C3(zFock(:,:,:,nspin),zSortedEigenvectors(:,:,nspin),ndim,RotateLayers,numb,numk,&
                DegFactor,nOccStates(nspin),nPartOccStates(nspin),Coords,MomentaValues,numNeighborCells,nUnitCell_1,nUnitCell_2,nSortedMomenta(:,:,nspin),nC3pairs)
            else
                call GetFock(zFock(:,:,:,nspin),zSortedEigenvectors(:,:,nspin),ndim,numb,numk,&
                DegFactor,nOccStates(nspin),nPartOccStates(nspin),Coords,MomentaValues,numNeighborCells,nUnitCell_1,nUnitCell_2,nSortedMomenta(:,:,nspin))
            endif

            call ValleyTransform(zEigenvectors(:,:,:,nspin),&
                Coords,MomentaValues,nMomentaComponents,numb,ndim,numk,Nk,numNeighborCells,nUnitCell_1,nUnitCell_2,reshape([t1,t2,t3],[2,3]),RotMatrix)
            
            call SortEnergies_1Spin(ndim,numk,numb,Nk,nMomentaComponents,zSortedEigenvectors(:,:,nspin),&
                SortedEnergies_1Spin(:,nspin),nSortedMomenta(:,:,nspin),zEigenvectors(:,:,:,nspin),Bands(:,:,nspin))

            if(nenforceC3.eq.1)then
                call GetFock_C3(zFock(:,:,:,nspin),zSortedEigenvectors(:,:,nspin),ndim,RotateLayers,numb,numk,&
                DegFactor,nOccStates(nspin),nPartOccStates(nspin),Coords,MomentaValues,numNeighborCells,nUnitCell_1,nUnitCell_2,nSortedMomenta(:,:,nspin),nC3pairs)
            else
                call GetFock(zFock(:,:,:,nspin),zSortedEigenvectors(:,:,nspin),ndim,numb,numk,&
                DegFactor,nOccStates(nspin),nPartOccStates(nspin),Coords,MomentaValues,numNeighborCells,nUnitCell_1,nUnitCell_2,nSortedMomenta(:,:,nspin))
            endif

            zFock(:,:,:,nspin)= .5_dp*zFock(:,:,:,nspin)+zFockBulk(:,:,:,nspin)
            zFock(:,:,:,nspin)=zFock(:,:,:,nspin)/cmplx(real(numk*numk,dp),0.0_dp,dp)
        
        else
            write(*,*) 'get Fock...'
            zFock(:,:,:,nspin) = cmplx(0.0_dp,0.0_dp,dp)
            if(nenforceC3.eq.1)then
                call GetFock_C3(zFock(:,:,:,nspin),zSortedEigenvectors(:,:,nspin),ndim,RotateLayers,numb,numk,&
                DegFactor,nOccStates(nspin),nPartOccStates(nspin),Coords,MomentaValues,numNeighborCells,nUnitCell_1,nUnitCell_2,nSortedMomenta(:,:,nspin),nC3pairs)
            else
                call GetFock(zFock(:,:,:,nspin),zSortedEigenvectors(:,:,nspin),ndim,numb,numk,&
                DegFactor,nOccStates(nspin),nPartOccStates(nspin),Coords,MomentaValues,numNeighborCells,nUnitCell_1,nUnitCell_2,nSortedMomenta(:,:,nspin))
            endif

            zFock(:,:,:,nspin) = zFock(:,:,:,nspin)+zFockBulk(:,:,:,nspin)
            zFock(:,:,:,nspin) = zFock(:,:,:,nspin)/cmplx(real(numk*numk,dp),0.0_dp,dp)
        endif

        ! Time reversal symmetry
        if(nenforceT.eq.1)then
            zFock(:,:,:,nspin) = cmplx(real(zFock(:,:,:,nspin),dp),0.0_dp,dp)
            zFockBulk(:,:,:,nspin) = cmplx(real(zFockBulk(:,:,:,nspin),dp),0.0_dp,dp)
        endif

        ! C2 symmetry
        if(nenforceC2.eq.1)then
            do m=1,numNeighborCells
                zFock(:,:,m,nspin) = .5_dp*(zFock(:,:,m,nspin) + conjg(transpose(zFock(nC2pairs(:),nC2pairs(:),m,nspin))))
                zFockBulk(:,:,m,nspin) = .5_dp*(zFockBulk(:,:,m,nspin) + conjg(transpose(zFockBulk(nC2pairs(:),nC2pairs(:),m,nspin))))
            enddo
        endif

        ! C2T symmetry
        if(nenforceC2T.eq.1)then
            do m=1,numNeighborCells
                zFock(:,:,m,nspin) = .5_dp*(zFock(:,:,m,nspin) + (transpose(zFock(nC2pairs(:),nC2pairs(:),m,nspin))))
                zFockBulk(:,:,m,nspin) = .5_dp*(zFockBulk(:,:,m,nspin) + (transpose(zFockBulk(nC2pairs(:),nC2pairs(:),m,nspin))))
            enddo
        endif

    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! Get Step, Fock, Potential, Energies !!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(*,*) 'Get Step, Fock, Potential, Energies'
    HartreeEnergyIn = HartreeEnergy
    FockEnergyIn = FockEnergy
    HubbardEnergyIn = HubbardEnergy
    KineticEnergyIn = KineticEnergy
    TotalEnergyIn = TotalEnergy

    do n=1,ndim
        Density(n,:) = real(zFock(n,n,1,:),dp)
    enddo

    call OptimalStep(Step,HartreeEnergy,FockEnergy,HubbardEnergy,KineticEnergy,zFock,zFockIn,&
        Density,DensityIn,DensitySub,ndim,numNeighborCells,Coords,LongRange,Delta,nUnitCell_1,nUnitCell_2,reshape([t1,t2,t3],[2,3]))

    TotalEnergy = alpha*sum(FockEnergy)*real(3-numS,dp) + sum(KineticEnergy)*real(3-numS,dp) + U*HubbardEnergy + alphaH*HartreeEnergy

    zFock = Step*zFock + (1.0_dp-Step)*zFockIn

    do n=1,ndim
        Density(n,:)=real(zFock(n,n,1,:),dp)
    enddo

    DensityConvergence(:) = 0.0_dp
    do n=1,ndim
        do nspin=1,numS
            if(DensityConvergence(nspin).LT.abs(Density(n,nspin)-DensityIn(n,nspin)))then
                DensityConvergence(nspin)=abs(Density(n,nspin)-DensityIn(n,nspin))
            endif
        enddo
    enddo

    DensityIn = Density
    zFockIn = zFock


    Potential=0.0_dp
    if(numS.EQ.1)then
        do n=1,ndim
        do m=1,ndim
            Potential(n,1) = Potential(n,1) + 2._dp*(Density(m,1)-DensitySub(m))*alphaH*LongRange(n,m)
        enddo
        enddo
        Potential(:,1) = Potential(:,1) + U*(Density(:,1)-DensitySub(:))
    endif
    if(numS.EQ.2)then
        do n=1,ndim
            do m=1,ndim
                Potential(n,1)=Potential(n,1)+(Density(m,1)+Density(m,2)-2.0_dp*DensitySub(m))*alphaH*LongRange(n,m)
            enddo
        enddo

        Potential(:,2)=Potential(:,1)+U*(Density(:,1)-DensitySub(:))
        Potential(:,1)=Potential(:,1)+U*(Density(:,2)-DensitySub(:))
    endif

    write(*,*) 'Output'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! Output bands !!! !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do nspin=1,numS
        if(nPrintOutputs.EQ.1)then

            write(filename,'(A7,A5,A6,I0,A210)') dir,'Bands','-nspin',nspin,parameters
            open(99,file=filename,status='replace')
            call PlotBands(Bands(:,:,nspin),FermiEnergy,nMomentaFlattened,numb,numk,Nk,g1,g12)
            close(99)

        endif
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! Output convergence and energies !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(nConvergenceFile+1,fmt='(F5.1,3X)', advance="no") real(it,dp)
    write(nConvergenceFile+1,fmt='(F9.6,3X)', advance="no") maxval([DensityConvergence(1),DensityConvergence(2)])
    write(nConvergenceFile+1,fmt='(F5.3,3X)', advance="no") Step
    if(numS.EQ.1)then
        write(nConvergenceFile+1,fmt='(F24.9,3X)') TotalEnergyIn - TotalEnergy
    else
        write(nConvergenceFile+1,fmt='(F24.9,3X)') TotalEnergyIn - TotalEnergy
    endif

    if(numS.EQ.1)then
        write(nEnergyFile+1,fmt='(F5.1,3X)', advance="no") real(it,dp)
        write(nEnergyFile+1,fmt='(F24.9,3X)', advance="no") TotalEnergyIn
        write(nEnergyFile+1,fmt='(F24.9,3X)', advance="no") 2.0_dp*KineticEnergyIn(1)
        write(nEnergyFile+1,fmt='(F24.9,3X)', advance="no") 2.0_dp*FockEnergyIn(1)
        write(nEnergyFile+1,fmt='(F24.9,3X)', advance="no") HartreeEnergyIn
        write(nEnergyFile+1,fmt='(F24.9,3X)') HubbardEnergyIn
    else
        write(nEnergyFile+1,fmt='(F5.1,3X)', advance="no") real(it,dp)
        write(nEnergyFile+1,fmt='(F24.9,3X)', advance="no") TotalEnergyIn
        write(nEnergyFile+1,fmt='(F24.9,3X)', advance="no") KineticEnergyIn(1)+KineticEnergyIn(2)
        write(nEnergyFile+1,fmt='(F24.9,3X)', advance="no") FockEnergyIn(1)+FockEnergyIn(2)
        write(nEnergyFile+1,fmt='(F24.9,3X)', advance="no") HartreeEnergyIn
        write(nEnergyFile+1,fmt='(F24.9,3X)') HubbardEnergyIn
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!! Output Fock !!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    write(*,*) 'writing Fock ...'

    if (nWriteFock.EQ.1) then

        do nspin=1,numS
            write(filename,'(A9,A10,I0,A210)') dirFock,'Fock-nspin',nspin,parameters
            open(91, file = filename, status = 'replace', form='UNFORMATTED', &
                ACCESS='direct', recl = 2*dp, iostat=my_iostat, iomsg=my_iomsg)
            if(my_iostat /= 0) then
                write(*,*) 'Write Fock failed with iostat = ', my_iostat, ' iomsg = '//trim(my_iomsg)
            endif
            
            rcc=0
                    do i=1,ndim
                        do j=1,i
                            rcc = rcc+1
                            write(91,rec=rcc, iostat=my_iostat, iomsg=my_iomsg) zFock(i,j,1,nspin)
                        enddo
                    enddo
            do m=2,numNeighborCells
                    do i=1,ndim
                        do j=1,ndim
                            rcc = rcc+1
                            write(91,rec=rcc, iostat=my_iostat, iomsg=my_iomsg) zFock(i,j,m,nspin)
                        enddo
                    enddo
            enddo
            close(91)
  
        enddo

        
    endif

enddo  !while loop

close(nEnergyFile+1)
close(nConvergenceFile+1)


end program FORGE

!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!cc
!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!cc

 
 

 
