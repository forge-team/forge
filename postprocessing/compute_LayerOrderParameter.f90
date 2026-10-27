! Interlayer, inter-sublattice, inter-valley order from the saved Fock density.
! Implements corrected OrderParameterLayer Eqs. (63)-(83); see layer_orderparameter.md.
! K follows OrderParameter.f90: K0=(-4*pi/3,0), Kprime=-K.
! The first two loop channels use the three-neighbor averages of Main_OrderParamsFull.
program compute_LayerOrderParameter
    use Setup, only: dp, nlayers, ntheta, RotateLayers, ndim, numI, numS, a0, a1, a2, pi, dir
    use Geometry, only: NumberNeighborCells, OrderNeighborCells
    use PostProcessingInput, only: BuildGeometry, InputParameters, ReadFock
    implicit none

    integer(dp) :: ind, nspin, i, sub, layer, first, last, loop
    integer(dp) :: center2, centerCell(2)
    integer(dp), parameter :: nLoopSamples(3)=[3_dp,3_dp,1_dp]
    integer(dp), allocatable :: ind_j(:), ind_l(:), centerMap(:), centerShift(:,:)
    integer(dp), allocatable :: rowSites(:,:,:,:), colSites(:,:,:,:), fockCells(:,:,:,:)
    logical, allocatable :: reverseCells(:,:,:,:), centerSeen(:)
    real(dp) :: t1(2), t2(2), t3(2), vq1(2), vq12(2), ang_mat(2,2), layerRotation(2,2)
    real(dp) :: layerVectors(2,2,2), deltaVec(2,3,2), valleyK(2,2), r1(2), r2(2)
    real(dp) :: geometryTolerance, cellDet
    real(dp), allocatable :: coord(:,:), coord_ref(:,:), outputCenters(:,:)
    complex(dp), allocatable :: zFock(:,:,:,:), chiPhase(:)
    complex(dp) :: Delta21(3), Delta12(3), localRho(4), summedRho(4), sublatticeRho(4,2)
    complex(dp), parameter :: zi=cmplx(0.0_dp,1.0_dp,dp)
    character(210) :: parameters
    character(512) :: filename, summaryFilename
    character(16) :: spinText
    character(12), parameter :: component(4) = [character(12) :: &
        'AK2_BKp1', 'BK2_AKp1', 'AK1_BKp2', 'BK1_AKp2']
    character(1) :: sublattice
    character(32) :: argument
    logical :: checkGeometry
    integer :: unit, summaryUnit, ios

    if(nlayers /= 2) error stop 'compute_LayerOrderParameter requires nlayers=2'
    if(ntheta < 1) error stop 'compute_LayerOrderParameter requires ntheta>=1'
    if(any(abs(RotateLayers) /= 1)) error stop 'RotateLayers must contain only -1 or +1'
    call get_command_argument(1,argument)
    if(len_trim(argument)>0 .and. trim(argument)/='--check-geometry') then
        error stop 'Usage: compute_LayerOrderParameter [--check-geometry]'
    endif
    checkGeometry=(trim(argument)=='--check-geometry')

    allocate(coord(ndim,3),coord_ref(ndim,3),outputCenters(ndim/2,2))
    allocate(centerMap(ndim),centerShift(2,ndim),centerSeen(ndim),chiPhase(ndim/2))
    allocate(rowSites(18,3,2,ndim/2),colSites(18,3,2,ndim/2), &
             fockCells(18,3,2,ndim/2),reverseCells(18,3,2,ndim/2))
    call NumberNeighborCells(numI,ind)
    allocate(ind_j(ind),ind_l(ind))
    call OrderNeighborCells(numI,ind,ind_j,ind_l)
    call BuildGeometry(coord,t1,t2,t3,vq1,vq12,ang_mat,coord_ref)
    cellDet=t1(1)*t2(2)-t1(2)*t2(1)
    if(abs(cellDet)<spacing(1.0_dp)) error stop 'Degenerate moire lattice'
    geometryTolerance=1.0e-5_dp*a0+128.0_dp*spacing(1.0_dp)*max(1.0_dp,maxval(abs(coord_ref)))

    ! The manuscript's primitive labels are (Setup a2,Setup a1). This reflection
    ! keeps its alpha=exp(2*pi*i/3) AND the legacy OrderParameter valley labels.
    ! d2 is the A-to-B basis displacement; all three deltas are nearest neighbors.
    do layer=1,2
        if(RotateLayers(layer)==-1) then
            layerRotation=ang_mat
        else
            layerRotation=transpose(ang_mat)
        endif
        layerVectors(:,1,layer)=matmul(layerRotation,a2)
        layerVectors(:,2,layer)=matmul(layerRotation,a1)
        deltaVec(:,2,layer)=(layerVectors(:,1,layer)+layerVectors(:,2,layer))/3.0_dp
        deltaVec(:,1,layer)=deltaVec(:,2,layer)-layerVectors(:,2,layer)
        deltaVec(:,3,layer)=deltaVec(:,2,layer)-layerVectors(:,1,layer)
        valleyK(:,layer)=matmul(layerRotation,[-4.0_dp*pi/3.0_dp,0.0_dp])
    enddo

    ! Pair physical centers once, then locate each bond on the exact reference
    ! lattice. Relaxation changes orbital positions, not the site labels/topology.
    do sub=1,2
        call SiteRange(1_dp,sub,first,last)
        call MatchCenters(sub,first,last,centerMap,centerShift)
        centerSeen=.false.
        sublattice='A'
        if(sub==2) sublattice='B'
        do i=first,last
            center2=centerMap(i)
            centerCell=centerShift(:,i)
            if(centerSeen(center2)) error stop 'Layer-center mapping is not one-to-one'
            centerSeen(center2)=.true.
            r1=coord_ref(i,1:2)
            r2=coord_ref(center2,1:2)+real(centerCell(1),dp)*t1+real(centerCell(2),dp)*t2
            chiPhase(i)=exp(zi*(dot_product(valleyK(:,1),r1)+dot_product(valleyK(:,2),r2)))
            outputCenters(i,:)=0.5_dp*(coord(i,1:2)+coord(center2,1:2) &
                +real(centerCell(1),dp)*t1+real(centerCell(2),dp)*t2)
            call CacheLoopSites(i,1_dp,2_dp,1_dp,sublattice,r2,r1)
            call CacheLoopSites(i,2_dp,1_dp,2_dp,sublattice,r1,r2)
        enddo
        call SiteRange(2_dp,sub,first,last)
        if(any(.not.centerSeen(first:last))) error stop 'Some layer-2 centers were not paired'
    enddo
    if(checkGeometry) then
        write(*,'(A)') 'Layer-order geometry check passed: paired centers, bond topology and Fock translations.'
        stop
    endif

    ! Same saved-state suffix and packed reader as Main_OrderParamsFull.
    allocate(zFock(ndim,ndim,ind,numS))
    call InputParameters(parameters)
    write(*,'(A,A)') 'Reading Fock state ',trim(parameters)
    call ReadFock(zFock,ind,parameters)

    do nspin=1,numS
        write(spinText,'(I0)') nspin
        filename=trim(dir)//'LayerOrderParameterLocal-numb'//intString(ndim)// &
            '-nspin'//trim(spinText)//trim(parameters)
        summaryFilename=trim(dir)//'LayerOrderParameter-numb'//intString(ndim)// &
            '-nspin'//trim(spinText)//trim(parameters)
        open(newunit=unit,file=trim(filename),status='replace',action='write',iostat=ios)
        if(ios/=0) then
            write(*,'(A,A)') 'Cannot open output file: ',trim(filename)
            error stop 1
        endif
        write(unit,'(A)') '# Interlayer inter-sublattice inter-valley envelope products f_bra^* f_ket (not 3*f_bra^* f_ket).'
        write(unit,'(A)') '# K0=(-4*pi/3,0), Kprime=-K, rotated with each layer; phase uses unrelaxed reference coordinates.'
        write(unit,'(A)') '# Physical center coordinates are in units a=0.246 nm. No spin multiplicity applied.'
        write(unit,'(A)') '# Legacy three-neighbor averaging in raw loop channels 1/2; central channel 3 is unaveraged.'
        write(unit,'(A)') '# Columns: i_layer1 i_layer2 sublattice x_center y_center Re/Im AK2_BKp1 BK2_AKp1 AK1_BKp2 BK1_AKp2'
        sublatticeRho=cmplx(0.0_dp,0.0_dp,dp)
        do sub=1,2
            call SiteRange(1_dp,sub,first,last)
            sublattice='A'
            if(sub==2) sublattice='B'
            do i=first,last
                call CachedLoopChannels(i,1_dp,Delta21)
                call CachedLoopChannels(i,2_dp,Delta12)
                call ProjectLayerChannels(Delta21,Delta12,sublattice,chiPhase(i),localRho)
                sublatticeRho(:,sub)=sublatticeRho(:,sub)+localRho
                write(unit,'(2(I0,1X),A1,10(1X,ES24.16))') i,centerMap(i),sublattice,outputCenters(i,:), &
                    (real(localRho(loop),dp),aimag(localRho(loop)),loop=1,4)
            enddo
        enddo
        close(unit)
        ! All centers are sampled, not every third Kekule center. Average the
        ! A/B reconstructions once: do not insert an additional factor of three.
        summedRho=0.5_dp*(sublatticeRho(:,1)+sublatticeRho(:,2))
        open(newunit=summaryUnit,file=trim(summaryFilename),status='replace',action='write',iostat=ios)
        if(ios/=0) then
            write(*,'(A,A)') 'Cannot open output file: ',trim(summaryFilename)
            error stop 1
        endif
        write(summaryUnit,'(A)') '# Interlayer inter-sublattice inter-valley density; one honeycomb unit cell counted once.'
        write(summaryUnit,'(A)') '# rho=0.5*sum_over_all_A_and_B_centers(local_product); rho_per_cell=rho/(ndim/4).'
        write(summaryUnit,'(A)') '# Columns: component Re(rho) Im(rho) Re(rho_per_cell) Im(rho_per_cell)'
        do loop=1,4
            write(summaryUnit,'(A,4(1X,ES24.16))') trim(component(loop)), &
                real(summedRho(loop),dp),aimag(summedRho(loop)), &
                real(summedRho(loop),dp)/real(ndim/4,dp),aimag(summedRho(loop))/real(ndim/4,dp)
        enddo
        close(summaryUnit)
        write(*,'(A,A)') 'Wrote ',trim(filename)
        write(*,'(A,A)') 'Wrote ',trim(summaryFilename)
    enddo

contains

    function intString(value) result(text)
        integer(dp), intent(in) :: value
        character(32) :: buffer
        character(:), allocatable :: text
        write(buffer,'(I0)') value
        text=trim(buffer)
    end function intString

    subroutine SiteRange(layer,sublatticeIndex,lo,hi)
        integer(dp), intent(in) :: layer,sublatticeIndex
        integer(dp), intent(out) :: lo,hi
        integer(dp) :: quarter,half
        quarter=ndim/4; half=ndim/2
        if(layer == 1) then
            if(sublatticeIndex == 1) then
                lo=1; hi=quarter
            else
                lo=quarter+1; hi=half
            endif
        else
            if(sublatticeIndex == 1) then
                lo=half+1; hi=3*quarter
            else
                lo=3*quarter+1; hi=ndim
            endif
        endif
    end subroutine SiteRange

    subroutine MatchCenters(sublatticeIndex,lo1,hi1,pairMap,pairShift)
        integer(dp), intent(in) :: sublatticeIndex,lo1,hi1
        integer(dp), intent(inout) :: pairMap(ndim),pairShift(2,ndim)
        integer(dp) :: lo2,hi2,n,i,j,ix,iy,i0,j0,j1,jj,matchedRow
        integer(dp), allocatable :: p(:),way(:),edgeShift(:,:,:)
        logical, allocatable :: used(:)
        real(dp), allocatable :: cost(:,:),u(:),v(:),minv(:)
        real(dp) :: dx,dy,dist2,deltaCost,cur

        call SiteRange(2_dp,sublatticeIndex,lo2,hi2)
        n=hi1-lo1+1
        if(hi2-lo2+1 /= n) error stop 'Layer sublattices have different atom counts'
        allocate(cost(n,n),edgeShift(2,n,n),p(0:n),way(0:n),used(0:n))
        allocate(u(0:n),v(0:n),minv(0:n))

        do i=1,n
            do j=1,n
                cost(i,j)=huge(1.0_dp)
                do ix=-1,1
                    do iy=-1,1
                        dx=coord(lo2+j-1,1)+real(ix,dp)*t1(1)+real(iy,dp)*t2(1)-coord(lo1+i-1,1)
                        dy=coord(lo2+j-1,2)+real(ix,dp)*t1(2)+real(iy,dp)*t2(2)-coord(lo1+i-1,2)
                        dist2=dx*dx+dy*dy
                        if(dist2 < cost(i,j)) then
                            cost(i,j)=dist2
                            edgeShift(:,i,j)=[ix,iy]
                        endif
                    enddo
                enddo
            enddo
        enddo

        ! Minimum-cost one-to-one pairing across the periodic moire cell.
        ! This assigns the relative displacement zeta_i without duplicating
        ! boundary centers when two independent nearest-neighbor choices tie.
        u=0.0_dp; v=0.0_dp; p=0; way=0
        do i=1,n
            p(0)=i; j0=0; minv=huge(1.0_dp); used=.false.; way=0
            do
                used(j0)=.true.
                i0=p(j0); deltaCost=huge(1.0_dp); j1=0
                do j=1,n
                    if(.not.used(j)) then
                        cur=cost(i0,j)-u(i0)-v(j)
                        if(cur < minv(j)) then
                            minv(j)=cur; way(j)=j0
                        endif
                        if(minv(j) < deltaCost) then
                            deltaCost=minv(j); j1=j
                        endif
                    endif
                enddo
                if(j1 == 0) error stop 'Failed to pair layer centers'
                do jj=0,n
                    if(used(jj)) then
                        u(p(jj))=u(p(jj))+deltaCost
                        v(jj)=v(jj)-deltaCost
                    else
                        minv(jj)=minv(jj)-deltaCost
                    endif
                enddo
                j0=j1
                if(p(j0) == 0) exit
            enddo
            do while(j0 /= 0)
                j1=way(j0)
                p(j0)=p(j1)
                j0=j1
            enddo
        enddo

        do j=1,n
            matchedRow=p(j)
            pairMap(lo1+matchedRow-1)=lo2+j-1
            pairShift(:,lo1+matchedRow-1)=edgeShift(:,matchedRow,j)
        enddo
    end subroutine MatchCenters

    subroutine LoopOffsets(centerSublattice,braLayer,ketLayer,loop,rowOffset,colOffset)
        character(1), intent(in) :: centerSublattice
        integer(dp), intent(in) :: braLayer,ketLayer,loop
        real(dp), intent(out) :: rowOffset(2,6),colOffset(2,6)
        real(dp) :: s,a1b(2),a2b(2),a1k(2),a2k(2),d1b(2),d2b(2),d3b(2),d1k(2),d2k(2),d3k(2)

        s=1.0_dp
        if(centerSublattice == 'B') s=-1.0_dp
        a1b=layerVectors(:,1,braLayer); a2b=layerVectors(:,2,braLayer)
        a1k=layerVectors(:,1,ketLayer); a2k=layerVectors(:,2,ketLayer)
        d1b=deltaVec(:,1,braLayer); d2b=deltaVec(:,2,braLayer); d3b=deltaVec(:,3,braLayer)
        d1k=deltaVec(:,1,ketLayer); d2k=deltaVec(:,2,ketLayer); d3k=deltaVec(:,3,ketLayer)

        select case(loop)
        case(1)
            rowOffset(:,1)= s*d1b;                 colOffset(:,1)=0.0_dp
            rowOffset(:,2)=-s*a2b;                 colOffset(:,2)= s*d1k
            rowOffset(:,3)=-s*a2b+s*d3b;           colOffset(:,3)=-s*a2k
            rowOffset(:,4)=-s*a1b;                 colOffset(:,4)=-s*a2k+s*d3k
            rowOffset(:,5)=-s*a1b+s*d2b;           colOffset(:,5)=-s*a1k
            rowOffset(:,6)=0.0_dp;                 colOffset(:,6)=-s*a1k+s*d2k
        case(2)
            rowOffset(:,1)= s*d2b;                 colOffset(:,1)=0.0_dp
            rowOffset(:,2)= s*a1b;                 colOffset(:,2)= s*d2k
            rowOffset(:,3)= s*a1b+s*d1b;           colOffset(:,3)= s*a1k
            rowOffset(:,4)= s*a1b-s*a2b;           colOffset(:,4)= s*a1k+s*d1k
            rowOffset(:,5)= s*d1b;                 colOffset(:,5)= s*a1k-s*a2k
            rowOffset(:,6)=0.0_dp;                 colOffset(:,6)= s*d1k
        case(3)
            rowOffset(:,1)= s*d3b;                 colOffset(:,1)=0.0_dp
            rowOffset(:,2)=-s*a1b+s*a2b;           colOffset(:,2)= s*d3k
            rowOffset(:,3)=-s*a1b+s*a2b+s*d2b;     colOffset(:,3)=-s*a1k+s*a2k
            rowOffset(:,4)= s*a2b;                 colOffset(:,4)=-s*a1k+s*a2k+s*d2k
            rowOffset(:,5)= s*a2b+s*d1b;           colOffset(:,5)= s*a2k
            rowOffset(:,6)=0.0_dp;                 colOffset(:,6)= s*a2k+s*d1k
        case default
            error stop 'Invalid loop channel'
        end select
    end subroutine LoopOffsets

    subroutine LoopAverageShift(centerSublattice,braLayer,ketLayer,loop,sample,rowShift,colShift)
        character(1), intent(in) :: centerSublattice
        integer(dp), intent(in) :: braLayer,ketLayer,loop,sample
        real(dp), intent(out) :: rowShift(2),colShift(2)
        integer(dp) :: n1,n2
        ! In FORGE's rotated Setup vectors the extra shifts are:
        ! loop 1: a1+a2, 2*a1-a2; loop 2: 2*a1-a2, a1-2*a2.
        ! layerVectors stores (Setup a2,Setup a1), hence the coefficients below.
        ! Each K.shift is an integer multiple of 2*pi, so chiPhase is unchanged.
        if(loop<1.or.loop>3) error stop 'Invalid averaged loop channel'
        if(sample<1.or.sample>nLoopSamples(loop)) error stop 'Invalid loop sample'
        n1=0; n2=0
        if(loop==1.and.sample==2) then
            n1=1; n2=1
        else if((loop==1.and.sample==3).or.(loop==2.and.sample==2)) then
            n1=-1; n2=2
        else if(loop==2.and.sample==3) then
            n1=-2; n2=1
        endif
        if(centerSublattice=='B') then
            n1=-n1; n2=-n2
        endif
        rowShift=real(n1,dp)*layerVectors(:,1,braLayer)+real(n2,dp)*layerVectors(:,2,braLayer)
        colShift=real(n1,dp)*layerVectors(:,1,ketLayer)+real(n2,dp)*layerVectors(:,2,ketLayer)
    end subroutine LoopAverageShift

    subroutine LocateReferenceSite(target,layer,sublatticeIndex,site,imageCell)
        real(dp), intent(in) :: target(2)
        integer(dp), intent(in) :: layer,sublatticeIndex
        integer(dp), intent(out) :: site,imageCell(2)
        integer(dp) :: lo,hi,j,candidateCell(2)
        real(dp) :: difference(2),fraction(2),residual(2),distance,bestDistance
        call SiteRange(layer,sublatticeIndex,lo,hi)
        bestDistance=huge(1.0_dp); site=0; imageCell=0
        do j=lo,hi
            difference=target-coord_ref(j,1:2)
            fraction(1)=(t2(2)*difference(1)-t2(1)*difference(2))/cellDet
            fraction(2)=(-t1(2)*difference(1)+t1(1)*difference(2))/cellDet
            candidateCell=nint(fraction,kind=dp)
            residual=difference-real(candidateCell(1),dp)*t1-real(candidateCell(2),dp)*t2
            distance=sqrt(dot_product(residual,residual))
            if(distance<bestDistance) then
                bestDistance=distance; site=j; imageCell=candidateCell
            endif
        enddo
        if(bestDistance>geometryTolerance) then
            write(*,*) 'No reference-lattice atom at loop endpoint; layer, sublattice, residual:', &
                layer,sublatticeIndex,bestDistance
            error stop 1
        endif
    end subroutine LocateReferenceSite

    subroutine CacheLoopSites(centerIndex,directionIndex,braLayer,ketLayer,centerSublattice,rBra,rKet)
        integer(dp), intent(in) :: centerIndex,directionIndex,braLayer,ketLayer
        character(1), intent(in) :: centerSublattice
        real(dp), intent(in) :: rBra(2),rKet(2)
        real(dp) :: rowOffset(2,6),colOffset(2,6),rowShift(2),colShift(2)
        integer(dp) :: m,n,sample,edge,centerSub,braSub,ketSub,rowSite,colSite,rowCell(2),colCell(2),shift(2),index
        logical :: reverse
        centerSub=1
        if(centerSublattice=='B') centerSub=2
        do m=1,3
            call LoopOffsets(centerSublattice,braLayer,ketLayer,m,rowOffset,colOffset)
            do sample=1,nLoopSamples(m)
            call LoopAverageShift(centerSublattice,braLayer,ketLayer,m,sample,rowShift,colShift)
            do n=1,6
                braSub=centerSub
                if(mod(n,2_dp)==1) braSub=3-centerSub
                ketSub=3-braSub
                call LocateReferenceSite(rBra+rowOffset(:,n)+rowShift,braLayer,braSub,rowSite,rowCell)
                call LocateReferenceSite(rKet+colOffset(:,n)+colShift,ketLayer,ketSub,colSite,colCell)
                shift=rowCell-colCell
                call FindFockCell(shift(1),shift(2),index)
                reverse=(index<1)
                if(reverse) call FindFockCell(-shift(1),-shift(2),index)
                if(index<1) then
                    write(*,*) 'Required Fock translation is outside Setup numI:',shift
                    error stop 1
                endif
                edge=6*(sample-1)+n
                rowSites(edge,m,directionIndex,centerIndex)=rowSite
                colSites(edge,m,directionIndex,centerIndex)=colSite
                fockCells(edge,m,directionIndex,centerIndex)=index
                reverseCells(edge,m,directionIndex,centerIndex)=reverse
            enddo
            enddo
        enddo
    end subroutine CacheLoopSites

    subroutine CachedLoopChannels(centerIndex,directionIndex,values)
        integer(dp), intent(in) :: centerIndex,directionIndex
        complex(dp), intent(out) :: values(3)
        integer(dp) :: m,n,rowSite,colSite,index
        complex(dp) :: overlap
        values=cmplx(0.0_dp,0.0_dp,dp)
        do m=1,3
            do n=1,6*nLoopSamples(m)
                rowSite=rowSites(n,m,directionIndex,centerIndex)
                colSite=colSites(n,m,directionIndex,centerIndex)
                index=fockCells(n,m,directionIndex,centerIndex)
                if(reverseCells(n,m,directionIndex,centerIndex)) then
                    overlap=conjg(zFock(colSite,rowSite,index,nspin))
                else
                    overlap=zFock(rowSite,colSite,index,nspin)
                endif
                values(m)=values(m)+overlap
            enddo
            values(m)=values(m)/real(nLoopSamples(m),dp)
        enddo
    end subroutine CachedLoopChannels

    subroutine ProjectLayerChannels(loops21,loops12,centerSublattice,chi,products)
        complex(dp), intent(in) :: loops21(3),loops12(3),chi
        character(1), intent(in) :: centerSublattice
        complex(dp), intent(out) :: products(4)
        complex(dp), parameter :: omega=cmplx(-0.5_dp,sqrt(3.0_dp)/2.0_dp,dp)
        complex(dp) :: ordered21(3),ordered12(3),a21,b21,a12,b12,prefactor
        ordered21=loops21; ordered12=loops12
        if(centerSublattice=='B') then
            ordered21=loops21([1,3,2]); ordered12=loops12([1,3,2])
        endif
        a21=(omega*ordered21(1)+ordered21(2)+conjg(omega)*ordered21(3))/3.0_dp
        b21=(conjg(omega)*ordered21(1)+ordered21(2)+omega*ordered21(3))/3.0_dp
        a12=(omega*ordered12(1)+ordered12(2)+conjg(omega)*ordered12(3))/3.0_dp
        b12=(conjg(omega)*ordered12(1)+ordered12(2)+omega*ordered12(3))/3.0_dp
        ! Eqs. (76)-(79) have 3*f_bra^* f_ket on the left. Return f_bra^* f_ket.
        ! The full sum phase is identical in the two layer directions.
        prefactor=zi*chi/(3.0_dp*sqrt(3.0_dp))
        products(1)= prefactor*(a21-omega*conjg(b12))
        products(2)=-prefactor*(omega*a21-conjg(b12))
        products(3)= prefactor*(a12-omega*conjg(b21))
        products(4)=-prefactor*(omega*a12-conjg(b21))
    end subroutine ProjectLayerChannels

    subroutine FindFockCell(r1,r2,index)
        integer(dp), intent(in) :: r1,r2
        integer(dp), intent(out) :: index
        integer(dp) :: k
        index=-1
        do k=1,ind
            if(ind_j(k) == r1 .and. ind_l(k) == r2) then
                index=k
                return
            endif
        enddo
    end subroutine FindFockCell


end program compute_LayerOrderParameter
