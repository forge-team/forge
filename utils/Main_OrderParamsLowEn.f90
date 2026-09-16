!
! Equivalent compute order parameters driver for the current FORGE workspace.
! This program reads an existing Fock matrix, constructs the real-space
! geometry, and computes order-parameter quantities using the external
! OrderPlot module routines.
!
program orderparams_lowenergy

use omp_lib
use Setup
use Geometry
use TightBinding
use HartreeFock
use OrderPlot
use lapack_routines

implicit none

character(300) :: filename
character(150) :: my_iomsg

integer(dp) :: Nk, icount, n, m, i, j, nspin, numNeighborCells, nband
integer(dp) :: nFermiLevel, nOccStates(numS), nPartOccStates(numS), nWindowSort_2Spins(2,2)
integer(dp) :: rcc, ivk1, ivk2

real(dp) :: vk(2)
real(dp), allocatable :: Energies(:)
complex(dp), allocatable :: zH(:,:)

real(dp) :: aMoire, cs, sn, FermiEnergy, DegFactor
real(dp) :: t1(2), t2(2), t3(2), tn(6,2), g1(2), g12(2), RotMatrix(2,2)
integer(dp) :: TnTonUnitCell12(0:6,2)

complex(dp) :: zinput

integer(dp), allocatable :: nMomentaComponents(:,:), nMomentaFlattened(:,:)
integer(dp), allocatable :: nSortedMomenta(:,:,:)
integer(dp), allocatable :: nUnitCell_1(:), nUnitCell_2(:)
integer(dp), allocatable :: NearestNeighborsUC(:,:), NearestNeighborsT(:,:)
integer(dp), allocatable :: nC2pairs(:), nC3pairs(:)
integer(dp), allocatable :: KekuleNeighbors(:,:,:), KekuleLattice(:,:,:)

real(dp), allocatable :: Coords(:,:)
real(dp), allocatable :: MomentaValues(:,:,:)
real(dp), allocatable :: Bands(:,:,:), SortedEnergies_1Spin(:,:), SortedEnergies_2Spins(:,:)
real(dp), allocatable :: Density(:,:), DensitySub(:)
real(dp), allocatable :: LongRange(:,:), Potential(:,:)
real(dp), allocatable :: ValleyPol(:), NormSquared(:)

complex(dp), allocatable :: zFockBulk(:,:,:,:), zFock(:,:,:,:)
complex(dp), allocatable :: zEigenvectors(:,:,:,:), zSortedEigenvectors(:,:,:)
complex(dp), allocatable :: fKA_conjxfKpB(:), fKB_conjxfKpA(:), fKp_conjxfK(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Sizes and allocations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Nk = numk*numk

call NumberNeighborCells(numI, numNeighborCells)

allocate(nMomentaComponents(Nk, 2))
allocate(nMomentaFlattened(numk, numk))
allocate(MomentaValues(1:numk, 1:numk, 1:2))
allocate(Coords(ndim, 3))
allocate(nC2pairs(ndim))
allocate(nC3pairs(ndim))
allocate(NearestNeighborsUC(ndim, 3))
allocate(NearestNeighborsT(ndim, 3))
allocate(Density(1:ndim, 1:numS))
allocate(DensitySub(ndim))
allocate(Potential(1:ndim, 1:numS))
allocate(LongRange(1:ndim, 1:ndim))

allocate(zEigenvectors(1:ndim, 1:numb, 1:Nk, 1:numS))
allocate(Bands(1:numb, 1:Nk, 1:numS))
allocate(zSortedEigenvectors(1:ndim, 1:numb*Nk, 1:numS))
allocate(nSortedMomenta(1:numb*Nk, 1:2, 1:numS))
allocate(SortedEnergies_1Spin(1:numb*Nk, 1:numS))
if(numS.EQ.2)then
    allocate(SortedEnergies_2Spins(1:2*numb*numk*numk,1:2))
endif

allocate(nUnitCell_1(1:numNeighborCells))
allocate(nUnitCell_2(1:numNeighborCells))

allocate(zFockBulk(1:ndim, 1:ndim, 1:numNeighborCells, 1:numS))
allocate(zFock(1:ndim, 1:ndim, 1:numNeighborCells, 1:numS))
zFockBulk(:,:,:,:) = cmplx(0.0_dp, 0.0_dp, dp)
zFock(:,:,:,:)     = cmplx(0.0_dp, 0.0_dp, dp)

Density(:,:)  = 0.0_dp
DensitySub(:) = 0.5_dp

allocate(KekuleNeighbors(ndim/2, 2, 3))
allocate(KekuleLattice(ndim/2, 6, 2))
allocate(fKA_conjxfKpB(ndim/2))
allocate(fKB_conjxfKpA(ndim/2))
allocate(fKp_conjxfK(ndim))
allocate(ValleyPol(ndim))
allocate(NormSquared(ndim))
allocate(Energies(ndim))
allocate(zH(ndim,ndim))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Geometry: reciprocal space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call OrderNeighborCells(numI, numNeighborCells, nUnitCell_1, nUnitCell_2)

aMoire = 3.0_dp*ntheta**2 + 3.0_dp*ntheta + 1.0_dp
g1  = (4.0_dp*pi/3.0_dp)/aMoire*(real(3*ntheta+1,dp)*a1+a2)
g12 = (4.0_dp*pi/3.0_dp)/aMoire*(real(3*ntheta+2,dp)*a2-a1)

cs = 1.0_dp-1.0_dp/(2.0_dp*aMoire)
sn = sqrt(1.0_dp-cs**2)
RotMatrix = reshape([cos(0.5_dp*acos(cs)),-sin(0.5_dp*acos(cs)),&
                     sin(0.5_dp*acos(cs)), cos(0.5_dp*acos(cs))],[2,2])

g1  = matmul(RotMatrix, g1)
g12 = matmul(RotMatrix, g12)

call SampleBZ(nMomentaComponents, nMomentaFlattened, MomentaValues, numk, g1, g12)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Geometry: real space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

t1 = real( ntheta  ,dp)*a1 + real(  ntheta+1,dp)*a2
t2 = real(-ntheta-1,dp)*a1 + real(2*ntheta+1,dp)*a2
t3 = t2-t1

call WignerSeitzCell(Coords, t1, t2, cs, sn)

do n=1,ndim
    Coords(n,:) = matmul(RotMatrix, Coords(n,1:2))
end do

t1 = matmul(RotMatrix, t1)
t2 = matmul(RotMatrix, t2)
t3 = t2-t1

if(nrelax.EQ.1)then
    if(nlayers.EQ.2) call LatticeRelaxationKoshino(Coords, g1, g12)
endif
if(nrelax.EQ.2)then
    if(nlayers.EQ.2) call LatticeRelaxationCarr(Coords, g1, g12)
endif

call C2_RelatedPoints(nC2pairs, Coords, ndim)
call C3_RelatedPoints(nC3pairs, Coords, ndim)

tn(1,:) = t1;  tn(2,:) = t2;  tn(3,:) = t3
tn(4,:) = -t1; tn(5,:) = -t2; tn(6,:) = -t3

TnTonUnitCell12(0,:) = [ 0,  0]
TnTonUnitCell12(1,:) = [ 1,  0]
TnTonUnitCell12(2,:) = [ 0,  1]
TnTonUnitCell12(3,:) = [-1,  1]
TnTonUnitCell12(4,:) = [-1,  0]
TnTonUnitCell12(5,:) = [ 0, -1]
TnTonUnitCell12(6,:) = [ 1, -1]

call getNearestNeighbors(NearestNeighborsUC, NearestNeighborsT, Coords, ndim, tn)

call getKekuleLattice(Coords, ndim, RotMatrix, a1, a2, tn, KekuleLattice)
call getKekuleNeighbors(KekuleNeighbors, Coords, ndim, a1, a2, RotMatrix, tn)

call LongRangeInteraction(LongRange, Coords, ndim, t1, t2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  File name strings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call InitParameters()

write(*,*) 'parameters',   parametersIn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Read Fock matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    do nspin=1,numS
        write(*,*) 'reading Fock...'
        if(numS.eq.1_dp)then
            write(filename,'(A9,A10,I0,A210)') dirFock,'Fock-nspin',1,parametersIn
        else
            write(filename,'(A9,A10,I0,A210)') dirFock,'Fock-nspin',nspin,parametersIn
        endif
        open(12, file=filename, form='unformatted', status='old', access='direct', recl=dp*2)
        rcc=0
        do i=1,ndim
            do j=1,i-1
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Density and Potential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Diagonalize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do nspin=1,numS
    write(*,*) 'diagonalizing nspin=',nspin
    !$omp parallel do &
    !$omp private(icount,ivk1,ivk2,vk,zH,Energies,nband) &
    !$omp shared(nspin,Coords,Nk,Bands,zEigenvectors,Potential,zFock,nMomentaComponents,MomentaValues) &
    !$omp shared(NearestNeighborsUC,NearestNeighborsT,numNeighborCells,nUnitCell_1,nUnitCell_2,t1,t2,t3)
    do icount=1,Nk
        ivk1=nMomentaComponents(icount,1)
        ivk2=nMomentaComponents(icount,2)
        vk(:)=MomentaValues(ivk1+1,ivk2+1,:)
        
        call HamiltonianHartreeFock(zH,Coords,Potential(:,nspin),alpha,Delta,zFock(:,:,:,nspin),&
            nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,vk,reshape([t1,t2,t3],[2,3]),&
            NearestNeighborsUC,NearestNeighborsT)
        
        call diagonalize(zH,Energies,'V',nLower,nUpper)
        
        do nband=1,nUpper-nLower+1
            Bands(nband,icount,nspin)=Energies(nband)
        enddo
        do nband=1,nUpper-nLower+1
            zEigenvectors(:,nband,icount,nspin)=zH(:,nband)
        enddo
    enddo
    !$omp end parallel do
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Sort energies and find Fermi level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

        if(nPartOccStates(1).gt.0)then
            DegFactor = 1.0_dp/real(nPartOccStates(1),dp)*real(nFermiLevel+nPartOccStates(1)-nOccStates(1),dp)
        else
            DegFactor = 0.0_dp
        endif
        write(*,*) 'Sorted energies',nspin, nOccStates(1), nPartOccStates(1), DegFactor

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Get Fock from eigenstates (numb bands)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do nspin=1,numS
    zFock(:,:,:,nspin) = cmplx(0.0_dp,0.0_dp,dp)
    if(nenforceC3.eq.1)then
        call GetFock_C3(zFock(:,:,:,nspin),zSortedEigenvectors(:,:,nspin),ndim,RotateLayers,numb,numk,&
        DegFactor,nOccStates(nspin),nPartOccStates(nspin),Coords,MomentaValues,numNeighborCells,&
        nUnitCell_1,nUnitCell_2,nSortedMomenta(:,:,nspin),nC3pairs)
    else
        call GetFock(zFock(:,:,:,nspin),zSortedEigenvectors(:,:,nspin),ndim,numb,numk,&
        DegFactor,nOccStates(nspin),nPartOccStates(nspin),Coords,MomentaValues,numNeighborCells,&
        nUnitCell_1,nUnitCell_2,nSortedMomenta(:,:,nspin))
    endif
    zFock(:,:,:,nspin) = zFock(:,:,:,nspin)/cmplx(real(numk*numk,dp),0.0_dp,dp)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Order parameters: numb-band Fock
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do nspin=1,numS

    call InterSubInterValAlt(fKA_conjxfKpB,fKB_conjxfKpA,KekuleLattice,KekuleNeighbors,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A21,I0,A6,I0,A210)') dir,'InterSubInterVal-numb',numb,'-nspin',nspin,parametersIn
    open(98,file=filename,status='replace')
    do i=1,ndim/2
        write(98,'(4(ES12.5,3X))') real(fKA_conjxfKpB(i)), aimag(fKA_conjxfKpB(i)), real(fKB_conjxfKpA(i)), aimag(fKB_conjxfKpA(i))
    enddo
    close(98)

    call IntraSubInterValAlt(fKp_conjxfK,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,NearestNeighborsUC,NearestNeighborsT,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A21,I0,A6,I0,A210)') dir,'IntraSubInterVal-numb',numb,'-nspin',nspin,parametersIn
    open(98,file=filename,status='replace')
    do i=1,ndim
        write(98,'(2(ES12.5,3X))') real(fKp_conjxfK(i)), aimag(fKp_conjxfK(i))
    enddo
    close(98)

    call InterSubIntraValAlt(fKA_conjxfKpB,fKB_conjxfKpA,KekuleLattice,KekuleNeighbors,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A21,I0,A6,I0,A210)') dir,'InterSubIntraVal-numb',numb,'-nspin',nspin,parametersIn
    open(98,file=filename,status='replace')
    do i=1,ndim/2
        write(98,'(2(ES12.5,3X))') real(fKA_conjxfKpB(i)), aimag(fKA_conjxfKpB(i))
    enddo
    do i=1,ndim/2
        write(98,'(2(ES12.5,3X))') real(fKB_conjxfKpA(i)), aimag(fKB_conjxfKpA(i))
    enddo
    close(98)

    call IntraSubIntraValAlt(ValleyPol,NormSquared,ndim,NearestNeighborsUC,NearestNeighborsT,numNeighborCells,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A21,I0,A6,I0,A210)') dir,'IntraSubIntraVal-numb',numb,'-nspin',nspin,parametersIn
    open(98,file=filename,status='replace')
    do i=1,ndim
        write(98,'(2(ES12.5,3X))') ValleyPol(i), NormSquared(i)
    enddo
    close(98)

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Background to subtract
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do nspin=1,numS
    do nband=1,numb/2-2
        do icount=1,Nk
            zSortedEigenvectors(:,(nband-1)*Nk+icount,nspin) = zEigenvectors(:,nband,icount,nspin)
            nSortedMomenta((nband-1)*Nk+icount,1,nspin) = nMomentaComponents(icount,1)
            nSortedMomenta((nband-1)*Nk+icount,2,nspin) = nMomentaComponents(icount,2)
        enddo
    enddo
enddo

do nspin=1,numS
    zFockBulk(:,:,:,nspin) = zFock(:,:,:,nspin)

    zFock(:,:,:,nspin) = cmplx(0.0_dp,0.0_dp,dp)
    if(nenforceC3.eq.1)then
        call GetFock_C3(zFock(:,:,:,nspin),zSortedEigenvectors(:,:,nspin),ndim,RotateLayers,numb,numk,&
        0.0_dp,(numb/2-2)*Nk,0_dp,Coords,MomentaValues,numNeighborCells,&
        nUnitCell_1,nUnitCell_2,nSortedMomenta(:,:,nspin),nC3pairs)
    else
        call GetFock(zFock(:,:,:,nspin),zSortedEigenvectors(:,:,nspin),ndim,numb,numk,&
        0.0_dp,(numb/2-2)*Nk,0_dp,Coords,MomentaValues,numNeighborCells,&
        nUnitCell_1,nUnitCell_2,nSortedMomenta(:,:,nspin))
    endif

    zFock(:,:,:,nspin) = zFockBulk(:,:,:,nspin) - zFock(:,:,:,nspin)/cmplx(real(numk*numk,dp),0.0_dp,dp)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Order parameters: 4-band contribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do nspin=1,numS

    call InterSubInterValAlt(fKA_conjxfKpB,fKB_conjxfKpA,KekuleLattice,KekuleNeighbors,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A28,I0,A210)') dir,'InterSubInterVal-numb4-nspin',nspin,parametersIn
    open(98,file=filename,status='replace')
    do i=1,ndim/2
        write(98,'(4(ES12.5,3X))') real(fKA_conjxfKpB(i)), aimag(fKA_conjxfKpB(i)), real(fKB_conjxfKpA(i)), aimag(fKB_conjxfKpA(i))
    enddo
    close(98)

    call IntraSubInterValAlt(fKp_conjxfK,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,NearestNeighborsUC,NearestNeighborsT,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A28,I0,A210)') dir,'IntraSubInterVal-numb4-nspin',nspin,parametersIn
    open(98,file=filename,status='replace')
    do i=1,ndim
        write(98,'(2(ES12.5,3X))') real(fKp_conjxfK(i)), aimag(fKp_conjxfK(i))
    enddo
    close(98)

    call InterSubIntraValAlt(fKA_conjxfKpB,fKB_conjxfKpA,KekuleLattice,KekuleNeighbors,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A28,I0,A210)') dir,'InterSubIntraVal-numb4-nspin',nspin,parametersIn
    open(98,file=filename,status='replace')
    do i=1,ndim/2
        write(98,'(2(ES12.5,3X))') real(fKA_conjxfKpB(i)), aimag(fKA_conjxfKpB(i))
    enddo
    do i=1,ndim/2
        write(98,'(2(ES12.5,3X))') real(fKB_conjxfKpA(i)), aimag(fKB_conjxfKpA(i))
    enddo
    close(98)

    call IntraSubIntraValAlt(ValleyPol,NormSquared,ndim,NearestNeighborsUC,NearestNeighborsT,numNeighborCells,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock(:,:,:,nspin))
    write(filename,'(A7,A28,I0,A210)') dir,'IntraSubIntraVal-numb4-nspin',nspin,parametersIn
    open(98,file=filename,status='replace')
    do i=1,ndim
        write(98,'(2(ES12.5,3X))') ValleyPol(i), NormSquared(i)
    enddo
    close(98)

enddo

end program orderparams_lowenergy

!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!cc
!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!c!cc
