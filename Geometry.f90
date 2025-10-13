module Geometry

use omp_lib
use Setup
use TightBinding

implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NumberNeighborCells(numI,numNeighborCells)

  integer(dp), intent(in) :: numI
  integer(dp), intent(inout) :: numNeighborCells
  integer(dp) :: n1, n2

  numNeighborCells=0_dp

  do n1=-numI,numI
    do n2=-numI,numI

        if(n1+n2.LE.numI.AND.n1+n2.GE.-numI)then
          if(n1.NE.0_dp.OR.n2.NE.0_dp)then
            numNeighborCells=numNeighborCells+1_dp
          endif
        endif

    enddo
  enddo
  numNeighborCells=numNeighborCells/2+1 ! only take half +1 due to symmetry Fock(i,j,n1,n2)=congj(Fock(j,i,-n1,-n2))(arg1,  arg2)

end subroutine NumberNeighborCells

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine OrderNeighborCells(numI, numNeighborCells, nUnitCell_1, nUnitCell_2)

  integer(dp), intent(in) :: numI, numNeighborCells
  integer(dp), intent(inout) :: nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells)

  integer(dp) :: n1, n2, icount
  
  nUnitCell_1=0_dp
  nUnitCell_2=0_dp
  icount=1

  nUnitCell_1(icount)=0_dp
  nUnitCell_2(icount)=0_dp

  if(numI.GT.0)then
    icount=icount+1
    nUnitCell_1(icount)=-1_dp
    nUnitCell_2(icount)=0_dp

    icount=icount+1
    nUnitCell_1(icount)=-1_dp
    nUnitCell_2(icount)=1_dp

    icount=icount+1
    nUnitCell_1(icount)=0_dp
    nUnitCell_2(icount)=-1_dp
  endif

  if(numI.GT.1)then
      do n1=-numI,numI
        do n2=-numI,numI
          if(n1+n2.LE.numI.AND.n1+n2.GE.-numI.AND.icount.LT.numNeighborCells)then
            if(n1.NE.0.OR.n2.NE.0)then
              icount=icount+1
              nUnitCell_1(icount)=n1
              nUnitCell_2(icount)=n2
              if(n1.EQ.-1.AND.n2.EQ.0)then
                  icount=icount-1
              endif
              if(n1.EQ.-1.AND.n2.EQ.1)then
                  icount=icount-1
              endif
              if(n1.EQ.0.AND.n2.EQ.-1)then
                  icount=icount-1
              endif
            endif
          endif
        enddo
      enddo
  endif

end subroutine OrderNeighborCells

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine SampleBZ(nMomentaComponents,nMomentaFlattened,nMomentaValues,numk,g1,g2)

  integer(dp), intent(in) :: numk
  real(dp), intent(in) :: g1(2), g2(2)
  integer(dp), intent(inout) :: nMomentaComponents(numk*numk,2),nMomentaFlattened(numk,numk)
  real(dp), intent(inout) :: nMomentaValues(numk,numk,2)

  integer(dp) :: ivk1,ivk2,icount

  icount=0
  do ivk1=0,numk-1
    do ivk2=0,numk-1
      icount=icount+1

      nMomentaComponents(icount,1)=ivk1
      nMomentaComponents(icount,2)=ivk2

      nMomentaFlattened(nMomentaComponents(icount,1)+1,nMomentaComponents(icount,2)+1)=icount
      
    enddo
  enddo

  do ivk2=0,numk-1
      do ivk1=0,numk-1
        nMomentaValues(ivk1+1,ivk2+1,1)=(ivk1*g1(1)+ivk2*g2(1))/real(numk,dp)
        nMomentaValues(ivk1+1,ivk2+1,2)=(ivk1*g1(2)+ivk2*g2(2))/real(numk,dp)
    enddo
  enddo
  


end subroutine SampleBZ

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine WignerSeitzCell(Coords,t1,t2,cs,sn)


  real(dp), intent(in)    :: t1(2), t2(2), cs, sn
  real(dp), intent(inout) :: Coords(ndim,3)
  integer(dp) :: n1, n2, n, nind, nrad, nlayer
  real(dp)    :: rMax, rTemp1, rTemp2, rTemp3, an1, an2, bn1, bn2, sq3=sqrt(3.0_dp)

  nrad=3*ntheta
  rMax = 3.0_dp*ntheta**2+3.0_dp*ntheta+1.0_dp

  nind=0
  do nlayer=1,nlayers

    if(RotateLayers(nlayer).eq.-1)then
      ! Find Coordinate of site A layer 1 is within 1st Wigner-Seitz cell
      do n1=-nrad,nrad
        do n2=-nrad,nrad
    
          an1=(n1+1.0_dp/3.0_dp)
          an2=(n2-2.0_dp/3.0_dp)
    
          rTemp1=abs(an1*(3.0_dp*ntheta+1)+an2*(3.0_dp*ntheta+2))+ 1e-6
          rTemp2=abs(an1-an2*(3.0_dp*ntheta+1))+ 1e-6
          rTemp3=abs(an1*(3.0_dp*ntheta+2)+an2)+ 1e-6

          if(rTemp1 < rMax .and. rTemp2 < rMax .and. rTemp3 < rMax)then

            nind=nind+1

            Coords(nind,1)=(n1+1.0_dp/3.0_dp)*a1(1)+(n2-2.0_dp/3.0_dp)*a2(1)
            Coords(nind,2)=(n1+1.0_dp/3.0_dp)*a1(2)+(n2-2.0_dp/3.0_dp)*a2(2)
            Coords(nind,3)= (real(nlayer-1) - real(nlayers-1)/2)*tz

          end if

        end do
      end do

      !! Include A atom at zone boundary
      nind=nind+1

      Coords(nind,1)=(-t2(1)-t1(1))/3.0_dp
      Coords(nind,2)=(-t2(2)-t1(2))/3.0_dp
      Coords(nind,3)= (real(nlayer-1) - real(nlayers-1)/2)*tz

      ! Find Coordsinate of site B layer 1 is within 1st Wigner-Seitz cell
      do n1=-nrad,nrad
        do n2=-nrad,nrad

          bn1=n1+2.0_dp/3.0_dp
          bn2=n2-1.0_dp/3.0_dp


          rTemp1=abs(bn1*(3.0_dp*ntheta+1)+bn2*(3.0_dp*ntheta+2))+ 1e-6
          rTemp2=abs(bn1-bn2*(3.0_dp*ntheta+1))+ 1e-6
          rTemp3=abs(bn1*(3.0_dp*ntheta+2)+bn2)+ 1e-6


          if(rTemp1 < rMax .and. rTemp2 < rMax .and. rTemp3 < rMax)then
            nind=nind+1

            Coords(nind,1)=bn1*a1(1)+bn2*a2(1)
            Coords(nind,2)=bn1*a1(2)+bn2*a2(2)
            Coords(nind,3)= (real(nlayer-1) - real(nlayers-1)/2)*tz
          end if


        end do
      end do

      !! Include B atom at zone boundary
      nind=nind+1

      Coords(nind,1)=(t2(1)+t1(1))/3.0_dp
      Coords(nind,2)=(t2(2)+t1(2))/3.0_dp
      Coords(nind,3)= (real(nlayer-1) - real(nlayers-1)/2)*tz

    else if (RotateLayers(nlayer).eq.1)then
      ! Find Coordsinate of site A layer 2 is within 1st Wigner-Seitz cell
      do n1=-nrad,nrad
        do n2=-nrad,nrad

          an1=(n1+1.0_dp/3.0_dp)*(cs-sn/sq3)-2*(n2-2.0_dp/3.0_dp)*sn/sq3
          an2=(n2-2.0_dp/3.0_dp)*(cs+sn/sq3)+2*(n1+1.0_dp/3.0_dp)*sn/sq3

          rTemp1=abs(an1*(3.0_dp*ntheta+1)+an2*(3.0_dp*ntheta+2))+ 1e-6
          rTemp2=abs(an1-an2*(3.0_dp*ntheta+1))+ 1e-6
          rTemp3=abs(an1*(3.0_dp*ntheta+2)+an2)+ 1e-6

          if(rTemp1 < rMax .and. rTemp2 .LT. rMax .and. rTemp3 < rMax)then
            nind=nind+1
            Coords(nind,1)=an1*a1(1)+an2*a2(1)
            Coords(nind,2)=an1*a1(2)+an2*a2(2)
            Coords(nind,3)= (real(nlayer-1) - real(nlayers-1)/2)*tz

          end if
        end do
      end do

      !! Include A atom at zone boundary
      nind=nind+1

      Coords(nind,1)=(t1(1)+t2(1))/3.0_dp
      Coords(nind,2)=(t1(2)+t2(2))/3.0_dp
      Coords(nind,3)= (real(nlayer-1) - real(nlayers-1)/2)*tz

      ! Find Coordsinate of site B layer 2 is within 1st Wigner-Seitz cell
      do n1=-nrad,nrad
        do n2=-nrad,nrad

          bn1=(n1+2.0_dp/3.0_dp)*(cs-sn/sq3)-2*(n2-1.0_dp/3.0_dp)*sn/sq3
          bn2=(n2-1.0_dp/3.0_dp)*(cs+sn/sq3)+2*(n1+2.0_dp/3.0_dp)*sn/sq3

          rTemp1=abs(bn1*(3.0_dp*ntheta+1)+bn2*(3.0_dp*ntheta+2))+ 1e-6
          rTemp2=abs(bn1-bn2*(3.0_dp*ntheta+1))+ 1e-6
          rTemp3=abs(bn1*(3.0_dp*ntheta+2)+bn2)+ 1e-6

          if(rTemp1 < rMax .AND. rTemp2 < rMax .AND. rTemp3 < rMax)then

            nind=nind+1

            Coords(nind,1)=bn1*a1(1)+bn2*a2(1)
            Coords(nind,2)=bn1*a1(2)+bn2*a2(2)
            Coords(nind,3)= (real(nlayer-1) - real(nlayers-1)/2)*tz

          end if
    
        end do
      end do

      !!! Include B atom at zone boundary
      nind=nind+1
    
      Coords(nind,1)=(-t1(1)-t2(1))/3.0_dp
      Coords(nind,2)=(-t1(2)-t2(2))/3.0_dp
      Coords(nind,3)= (real(nlayer-1) - real(nlayers-1)/2)*tz
    
    endif
  enddo

end subroutine WignerSeitzCell

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine C2_RelatedPoints(nC2pairs,Coords,ndim)

  integer(dp), intent(in) :: ndim
  real(dp), intent(in) :: Coords(ndim,3)
  integer(dp), intent(inout) :: nC2pairs(ndim)

  integer(dp) :: i,j,cnt,ilayer
  real(dp) :: r

  cnt=0
  do ilayer=1,nlayers
    do i=ndim/nlayers*(ilayer-1) + 1, ndim/nlayers*(ilayer-1) + ndim/nlayers/2
      do j=ndim/nlayers*(ilayer-1) + ndim/nlayers/2 + 1, ndim/nlayers*(ilayer-1) + ndim/nlayers 
        r = sqrt((Coords(i,1)+Coords(j,1))**2+(Coords(i,2)+Coords(j,2))**2)
        if (r.lt.1e-5)then
          nC2pairs(i) = j
          nC2pairs(j) = i
          cnt = cnt+1
        endif
      enddo
    enddo
  enddo
  
  if(cnt.lt.ndim/2)then
    write(*,*) 'ERROR C2. cnt=',cnt
  endif

end subroutine C2_RelatedPoints 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine C3_RelatedPoints(nC3pairs,Coords,ndim)

  integer(dp), intent(in) :: ndim
  real(dp), intent(in) :: Coords(ndim,3)
  integer(dp), intent(inout) :: nC3pairs(ndim)

  integer(dp) :: i,j,cnt,nlayer
  real(dp) :: r,x,y,x3,y3

  cnt=0
  do nlayer=1,nlayers
    do i=(nlayer-1)*ndim/nlayers + 1, (nlayer-1)*ndim/nlayers + ndim/nlayers/2
      do j=(nlayer-1)*ndim/nlayers + 1, (nlayer-1)*ndim/nlayers + ndim/nlayers/2
        x = Coords(i,1)
        y = Coords(i,2)
        x3 = -.5_dp*x - .5_dp*sqrt(3.0_dp)*y
        y3 = .5_dp*sqrt(3.0_dp)*x - .5_dp*y
        r = sqrt((Coords(j,1)-x3)**2+(Coords(j,2)-y3)**2)
        if (r.lt.1e-3)then
          nC3pairs(i) = j
          cnt = cnt+1
        endif
      enddo
    enddo

    do i=(nlayer-1)*ndim/nlayers + ndim/nlayers/2 + 1, (nlayer-1)*ndim/nlayers + ndim/nlayers
      do j=(nlayer-1)*ndim/nlayers + ndim/nlayers/2 + 1, (nlayer-1)*ndim/nlayers + ndim/nlayers
        x = Coords(i,1)
        y = Coords(i,2)
        x3 = -.5_dp*x - .5_dp*sqrt(3.0_dp)*y
        y3 = .5_dp*sqrt(3.0_dp)*x - .5_dp*y
        r = sqrt((Coords(j,1)-x3)**2+(Coords(j,2)-y3)**2)
        if (r.lt.1e-3)then
          nC3pairs(i) = j
          cnt = cnt+1
        endif
      enddo
    enddo

  enddo

  do nlayer=1,nlayers
    nC3pairs((nlayer-1)*ndim/nlayers + ndim/nlayers/2) = (nlayer-1)*ndim/nlayers + ndim/nlayers/2
    nC3pairs((nlayer-1)*ndim/nlayers + ndim/nlayers) = (nlayer-1)*ndim/nlayers + ndim/nlayers
    cnt = cnt+2
  enddo

  if(cnt.lt.ndim)then
          write(*,*) 'ERROR C3. cnt=',cnt
  endif

end subroutine C3_RelatedPoints 


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine Mz_RelatedPoints(nMzpairs,Coords,ndim)

  integer(dp), intent(in) :: ndim
  real(dp), intent(in) :: Coords(ndim,3)
  integer(dp), intent(inout) :: nMzpairs(ndim)

  integer(dp) :: i,j,cnt,nlayer
  real(dp) :: r

  cnt=0
  if(modulo(nlayers,2).eq.0)then
    do nlayer=1,nlayers/2
      do i = 1, ndim/nlayers
        nMzpairs(ndim/nlayers*(nlayer-1) + i) = ndim/nlayers*(nlayers - nlayer) + i
        nMzpairs(ndim/nlayers*(nlayers - nlayer) + i) = ndim/nlayers*(nlayer-1) + i
        cnt = cnt +2
      enddo
    enddo
  else
    do nlayer=1,(nlayers-1)/2
      do i = 1, ndim/nlayers
        nMzpairs(ndim/nlayers*(nlayer-1) + i) = ndim/nlayers*(nlayers - nlayer) + i
        nMzpairs(ndim/nlayers*(nlayers - nlayer) + i) = ndim/nlayers*(nlayer-1) + i
        cnt = cnt +2
      enddo
    enddo

    nlayer = (nlayers+1)/2
    do i = 1, ndim/nlayers
      nMzpairs(ndim/nlayers*(nlayer-1) + i) = ndim/nlayers*(nlayer-1) + i
      cnt = cnt +1
    enddo

  endif

  if(cnt.lt.ndim)then
    write(*,*) 'ERROR Mz. cnt=',cnt
  endif

end subroutine Mz_RelatedPoints 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Geometry
