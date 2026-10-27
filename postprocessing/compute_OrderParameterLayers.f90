! Real-space interlayer loop overlaps from OrderParameterLayer.pdf, Eqs. (63)-(65).
! The later continuum valley projections in that note are not applied here: the
! printed conjugations in Eqs. (27)-(52) and (67) are inconsistent with the
! defining overlaps and need correction before those projected rho values are
! safe to report.

program compute_OrderParameterLayers
    use Setup
    use Geometry, only: NumberNeighborCells, OrderNeighborCells
    use PostProcessingInput, only: BuildGeometry, InputParameters, ReadFock
    implicit none

    integer(dp) :: numNeighborCells, nspin, i, first, last, sub
    integer(dp), allocatable :: nUnitCell_1(:), nUnitCell_2(:)
    real(dp) :: t1(2), t2(2), t3(2), g1(2), g12(2), RotMatrix(2,2)
    real(dp) :: cellVectors(0:6,2), layerVectors(2,2,2), deltaVec(2,3,2)
    real(dp), allocatable :: Coords(:,:)
    logical, allocatable :: centerSeen(:)
    integer(dp), allocatable :: centerMap(:), centerShift(:,:)
    complex(dp), allocatable :: zFock(:,:,:,:)
    complex(dp) :: Delta21(3), Delta12(3)
    integer(dp) :: center2, centerCell(2)
    real(dp) :: r1(2), r2(2), center(2)
    character(210) :: inputSuffix
    character(512) :: filename
    character(16) :: spinText
    character(1) :: sublattice
    character(32) :: argument
    logical :: checkGeometry
    integer :: unit, ios

    allocate(Coords(ndim,3))
    allocate(centerSeen(ndim))
    allocate(centerMap(ndim),centerShift(2,ndim))
    call NumberNeighborCells(numI,numNeighborCells)
    allocate(nUnitCell_1(numNeighborCells),nUnitCell_2(numNeighborCells))
    call OrderNeighborCells(numI,numNeighborCells,nUnitCell_1,nUnitCell_2)
    allocate(zFock(ndim,ndim,numNeighborCells,numS))

    ! Keep geometry construction identical to the migrated band/order-parameter drivers.
    call BuildGeometry(Coords,t1,t2,t3,g1,g12,RotMatrix)
    cellVectors(0,:) = [0.0_dp,0.0_dp]
    cellVectors(1,:) = t1;  cellVectors(2,:) = t2;  cellVectors(3,:) = t3
    cellVectors(4,:) = -t1; cellVectors(5,:) = -t2; cellVectors(6,:) = -t3

    ! WignerSeitzCell rotates layer 1 by RotMatrix and layer 2 by its transpose.
    layerVectors(:,1,1) = matmul(RotMatrix,a1)
    layerVectors(:,2,1) = matmul(RotMatrix,a2)
    layerVectors(:,1,2) = matmul(transpose(RotMatrix),a1)
    layerVectors(:,2,2) = matmul(transpose(RotMatrix),a2)

    ! The PDF leaves delta's directions implicit in missing figures. This is
    ! the standard honeycomb basis convention used by Setup's a1,a2 ordering:
    ! d1=(a1+a2)/3, d2=d1+a1, d3=d1+a2. These are the three B positions
    ! used around a hexagon; they are not three independent bond vectors.
    do sub=1,2
        deltaVec(:,1,sub) = (layerVectors(:,1,sub)+layerVectors(:,2,sub))/3.0_dp
        deltaVec(:,2,sub) = deltaVec(:,1,sub)+layerVectors(:,1,sub)
        deltaVec(:,3,sub) = deltaVec(:,1,sub)+layerVectors(:,2,sub)
    enddo

    do sub=1,2
        call SiteRange(1_dp,sub,first,last)
        call MatchCenters(sub,first,last,centerMap,centerShift)
    enddo

    call InputParameters(inputSuffix)
    call get_command_argument(1,argument)
    checkGeometry=(trim(argument) == '--check-geometry')
    if(checkGeometry) then
        zFock=cmplx(0.0_dp,0.0_dp,dp)
        write(*,'(A)') 'Checking layer pairing, loop-site lookup, and required cell translations.'
    else
        write(*,'(A,A)') 'Reading Fock state ',trim(inputSuffix)
        call ReadFock(zFock,numNeighborCells,inputSuffix)
    endif

    do nspin=1,numS
        write(spinText,'(I0)') nspin
        filename = trim(dir)//'OrderParameterLayers-numb'//trim(intString(ndim))// &
            '-nspin'//trim(spinText)//trim(inputSuffix)
        if(.not.checkGeometry) then
            open(newunit=unit,file=trim(filename),status='replace',action='write',iostat=ios)
            if(ios /= 0) then
                write(*,'(A,A)') 'Cannot open output file: ',trim(filename)
                error stop 1
            endif
            write(unit,'(A)') '# Raw six-term interlayer overlap sums, Eqs. (63)-(65); no valley projection applied.'
            write(unit,'(A)') '# Assumed honeycomb convention: d1=(a1+a2)/3, d2=d1+a1, d3=d1+a2.'
            write(unit,'(A)') '# Fock convention: zFock(i,j,R) = sum_k psi_i(k)^* psi_j(k) exp(-i k.R).'
            write(unit,'(A)') '# Columns: center_index sublattice x_center y_center; Re/Im Delta21[m=1,2,3]; Re/Im Delta12[m=1,2,3].'
        endif

        do sub=1,2
            call SiteRange(1_dp,sub,first,last)
            centerSeen=.false.
            if(sub == 1) then
                sublattice='A'
            else
                sublattice='B'
            endif
            do i=first,last
                center2=centerMap(i)
                centerCell=centerShift(:,i)
                if(centerSeen(center2)) then
                    write(*,*) 'Layer-center mapping is not one-to-one:',center2
                    error stop 1
                endif
                centerSeen(center2)=.true.
                r1=Coords(i,1:2)
                r2=Coords(center2,1:2)+real(centerCell(1),dp)*t1+real(centerCell(2),dp)*t2
                center=0.5_dp*(r1+r2)

                call LoopChannels(2_dp,1_dp,sublattice,r2,r1,Delta21)
                call LoopChannels(1_dp,2_dp,sublattice,r1,r2,Delta12)
                if(.not.checkGeometry) then
                    write(unit,'(I0,1X,A1,14(1X,ES22.14))') i,sublattice,center, &
                        real(Delta21(1)),aimag(Delta21(1)),real(Delta21(2)),aimag(Delta21(2)), &
                        real(Delta21(3)),aimag(Delta21(3)),real(Delta12(1)),aimag(Delta12(1)), &
                        real(Delta12(2)),aimag(Delta12(2)),real(Delta12(3)),aimag(Delta12(3))
                endif
            enddo
            call SiteRange(2_dp,sub,first,last)
            if(any(.not.centerSeen(first:last))) then
                write(*,*) 'Some layer-2 centers were not paired for sublattice ',sublattice
                error stop 1
            endif
        enddo
        if(.not.checkGeometry) then
            close(unit)
            write(*,'(A,A)') 'Wrote ',trim(filename)
        endif
    enddo
    if(checkGeometry) write(*,'(A)') 'Geometry lookup check passed; no Fock state was read and no output was written.'

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
                        dx=Coords(lo2+j-1,1)+real(ix,dp)*t1(1)+real(iy,dp)*t2(1)-Coords(lo1+i-1,1)
                        dy=Coords(lo2+j-1,2)+real(ix,dp)*t1(2)+real(iy,dp)*t2(2)-Coords(lo1+i-1,2)
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

    subroutine LocateSite(target,lo,hi,site,imageCell,residual)
        real(dp), intent(in) :: target(2)
        integer(dp), intent(in) :: lo,hi
        integer(dp), intent(out) :: site,imageCell(2)
        real(dp), intent(out) :: residual
        integer(dp) :: j,q
        real(dp) :: d(2),dist

        residual=huge(1.0_dp); site=0; imageCell=0
        do j=lo,hi
            do q=0,6
                d=Coords(j,1:2)+cellVectors(q,:)-target
                dist=sqrt(dot_product(d,d))
                if(dist < residual) then
                    residual=dist; site=j
                    imageCell=[0,0]
                    if(q > 0) imageCell=CellCoordinates(q)
                endif
            enddo
        enddo
    end subroutine LocateSite

    function CellCoordinates(q) result(cell)
        integer(dp), intent(in) :: q
        integer(dp) :: cell(2)
        select case(q)
        case(1); cell=[ 1, 0]
        case(2); cell=[ 0, 1]
        case(3); cell=[-1, 1]
        case(4); cell=[-1, 0]
        case(5); cell=[ 0,-1]
        case(6); cell=[ 1,-1]
        case default; cell=[0,0]
        end select
    end function CellCoordinates

    subroutine LoopChannels(braLayer,ketLayer,centerSublattice,rBra,rKet,values)
        integer(dp), intent(in) :: braLayer,ketLayer
        character(1), intent(in) :: centerSublattice
        real(dp), intent(in) :: rBra(2),rKet(2)
        complex(dp), intent(out) :: values(3)
        real(dp) :: rowOffset(2,6),colOffset(2,6),rowTarget(2),colTarget(2),residual
        integer(dp) :: loop,term,rowCellLocal(2),colCellLocal(2)
        integer(dp) :: rowLocal,colLocal

        ! The lookup searches both sublattices because each endpoint is
        ! specified by its real-space position in the loop definition.
        values=cmplx(0.0_dp,0.0_dp,dp)
        do loop=1,3
            call LoopOffsets(centerSublattice,braLayer,ketLayer,loop,rowOffset,colOffset)
            do term=1,6
                rowTarget=rBra+rowOffset(:,term)
                colTarget=rKet+colOffset(:,term)
                call LocateLayerSite(rowTarget,braLayer,rowLocal,rowCellLocal,residual)
                call LocateLayerSite(colTarget,ketLayer,colLocal,colCellLocal,residual)
                values(loop)=values(loop)+FockElement(rowLocal,colLocal, &
                    rowCellLocal(1)-colCellLocal(1),rowCellLocal(2)-colCellLocal(2))
            enddo
        enddo
    end subroutine LoopChannels

    subroutine LocateLayerSite(target,layer,site,imageCell,residual)
        real(dp), intent(in) :: target(2)
        integer(dp), intent(in) :: layer
        integer(dp), intent(out) :: site,imageCell(2)
        real(dp), intent(out) :: residual
        integer(dp) :: loA,hiA,loB,hiB,siteA,siteB,cellA(2),cellB(2)
        real(dp) :: errorA,errorB
        call SiteRange(layer,1_dp,loA,hiA)
        call SiteRange(layer,2_dp,loB,hiB)
        call LocateSite(target,loA,hiA,siteA,cellA,errorA)
        call LocateSite(target,loB,hiB,siteB,cellB,errorB)
        if(errorA <= errorB) then
            site=siteA; imageCell=cellA; residual=errorA
        else
            site=siteB; imageCell=cellB; residual=errorB
        endif
        if(residual > 0.25_dp*a0) then
            write(*,*) 'No atom found at requested lattice position; residual, layer:',residual,layer
            error stop 1
        endif
    end subroutine LocateLayerSite

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

    function FockElement(rowSiteLocal,colSiteLocal,r1,r2) result(value)
        integer(dp), intent(in) :: rowSiteLocal,colSiteLocal,r1,r2
        complex(dp) :: value
        integer(dp) :: index

        call FindFockCell(r1,r2,index)
        if(index > 0) then
            value=zFock(rowSiteLocal,colSiteLocal,index,nspin)
            return
        endif
        call FindFockCell(-r1,-r2,index)
        if(index <= 0) then
            write(*,*) 'Required Fock translation is outside Setup numI:',r1,r2
            error stop 1
        endif
        ! The Fock files store one representative from each +/- cell pair.
        value=conjg(zFock(colSiteLocal,rowSiteLocal,index,nspin))
    end function FockElement

    subroutine FindFockCell(r1,r2,index)
        integer(dp), intent(in) :: r1,r2
        integer(dp), intent(out) :: index
        integer(dp) :: k
        index=-1
        do k=1,numNeighborCells
            if(nUnitCell_1(k) == r1 .and. nUnitCell_2(k) == r2) then
                index=k
                return
            endif
        enddo
    end subroutine FindFockCell

end program compute_OrderParameterLayers
