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

!!!!!!!!!!!!!!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine getNearestNeighbors(NearestNeighborsUC,NearestNeighborsT,Coords,ndim,tn)

    integer(dp), intent(in)    :: ndim
    integer(dp), intent(inout) :: NearestNeighborsUC(ndim,3),NearestNeighborsT(ndim,3)
    real(dp), intent(in)    :: tn(6,2), Coords(ndim,3)

    ! NearestNeighborsUC(i,j) = atom index inside the unit cell of the j^th nearest neighbor of the i^th atom
    ! NearestNeighborsT(i,j) = lattice vector displacement of the j^th nearest neighbor of the i^th atom (the
    !                           nearest neighbor lives outside the first unit cell fro atoms near the boundaries)
        
    integer(dp) :: i,j,n,ncount,ntemp,ntemp1,ntemp2,ntempt,ntempt1,ntempt2
    real(dp) :: rij,x,y,phi,temp,temp1,temp2,tempt,tempt1,tempt2,rot(ndim,3) 
    
    
    do i=1,ndim
        ncount=0
        ! write(*,*) i, Coords(i,1), Coords(i,2), Coords(i,3)
        do j=1,ndim
    
            if((Coords(i,3)-Coords(j,3))**2 .LT.tz*.01_dp)then
                
    
                rij=(Coords(i,1)-Coords(j,1))**2+(Coords(i,2)-Coords(j,2))**2
    
                rij=sqrt(rij)
    
                if(rij.GT.(a0*.2_dp).AND.rij.LT.(a0*1.2))then
                    ncount=ncount+1
                    NearestNeighborsUC(i,ncount)=j
                    NearestNeighborsT(i,ncount)=0
    
                    x=Coords(i,1)-Coords(j,1)
                    y=Coords(i,2)-Coords(j,2)
    
                    if(x.GT.0)then
                        phi=atan(y/x)
                    elseif(x.LT.0.AND.y.GT.0)then
                        phi=atan(y/x)+pi
                    else
                        phi=-pi+atan(y/x)   
                    endif
                    if(x.EQ.0)then
                        if(y.GT.0)then
                            phi=pi/2.
                        else
                            phi=-pi/2.
                        endif
                    endif
                    ! write(*,*) x,y,phi*180/pi
    
                    rot(i,ncount)=phi
                    ! write(*,*) i,ncount
                    ! write(*,*) Coords(i,3), Coords(j,3)
                    ! write(*,*) rij
                endif
    
                !!!!!!! Coupling to adjacent Wigner-Seitz cell
    
                do n=1,6
    
                    rij=(Coords(i,1)-Coords(j,1)+tn(n,1))**2+(Coords(i,2)-Coords(j,2)+tn(n,2))**2
    
                    rij=sqrt(rij)
    
                    if(rij.LT.a0*1.2)then
                        ncount=ncount+1
                        NearestNeighborsUC(i,ncount)=j
                        NearestNeighborsT(i,ncount)=n
    
                        x=Coords(i,1)-Coords(j,1)+tn(n,1)
                        y=Coords(i,2)-Coords(j,2)+tn(n,2)

                        if(x.GT.0)then
                            phi=atan(y/x)
                        elseif(x.LT.0.AND.y.GT.0)then
                            phi=atan(y/x)+pi
                        else
                            phi=-pi+atan(y/x)   
                        endif
                        if(x.EQ.0)then
                            if(y.GT.0)then
                                phi=pi/2.
                            else
                                phi=-pi/2.
                            endif
                        endif
    
    
                        rot(i,ncount)=phi
                        ! write(*,*) i,ncount
                        ! write(*,*) Coords(i,3), Coords(j,3)
                        ! write(*,*) rij
                    endif
    
                enddo
    
            endif
    
        enddo
    
        if(ncount.NE.3)then
            write(*,*) 'ERROR getNearestNeighbors', i,ncount
        endif
    
        !!!!!!!!!!!!!!!!!!!
        !!!!!! Order cycle
        !!!!!!!!!!!!!!!!!!!
    
    
        if(rot(i,1).LT.rot(i,2))then
            temp=rot(i,1)
            ntemp=NearestNeighborsUC(i,1)
            ntempt=NearestNeighborsT(i,1)
            rot(i,1)=rot(i,2)
            NearestNeighborsUC(i,1)=NearestNeighborsUC(i,2)
            NearestNeighborsT(i,1)=NearestNeighborsT(i,2)
            rot(i,2)=temp
            NearestNeighborsUC(i,2)=ntemp
            NearestNeighborsT(i,2)=ntempt
        endif
    
        if(rot(i,1).LT.rot(i,3))then
            temp1=rot(i,1)
            ntemp1=NearestNeighborsUC(i,1)
            ntempt1=NearestNeighborsT(i,1)
            temp2=rot(i,2)
            ntemp2=NearestNeighborsUC(i,2)
            ntempt2=NearestNeighborsT(i,2)
    
            rot(i,1)=rot(i,3)
            NearestNeighborsUC(i,1)=NearestNeighborsUC(i,3)
            NearestNeighborsT(i,1)=NearestNeighborsT(i,3)
            rot(i,2)=temp1
            NearestNeighborsUC(i,2)=ntemp1
            NearestNeighborsT(i,2)=ntempt1
            rot(i,3)=temp2
            NearestNeighborsUC(i,3)=ntemp2
            NearestNeighborsT(i,3)=ntempt2
            else
            if(rot(i,2).LT.rot(i,3))then
                temp=rot(i,2)
                ntemp=NearestNeighborsUC(i,2)
                ntempt=NearestNeighborsT(i,2)
                rot(i,2)=rot(i,3)
                NearestNeighborsUC(i,2)=NearestNeighborsUC(i,3)
                NearestNeighborsT(i,2)=NearestNeighborsT(i,3)
                rot(i,3)=temp
                NearestNeighborsUC(i,3)=ntemp
                NearestNeighborsT(i,3)=ntempt
            endif
        endif
    
    
    enddo
     
    
    return
end subroutine getNearestNeighbors


! !!!!!!!!!!!!!!!
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
      ! Find coordinate of site A layer 1 is within 1st Wigner-Seitz cell
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

      ! Find coordinate of site B layer 2 is within 1st Wigner-Seitz cell
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


subroutine LatticeRelaxationKoshino(Coords,g1,g12)

  real(dp), intent(inout) :: Coords(ndim,3)
  real(dp), intent(in)    :: g1(2), g12(2) 

  real(dp) :: Gn(8,2)
  integer(dp) :: i,j, m, n, Nk, h
  real(dp) :: Coordsnr(ndim,3), Coordstemp(ndim,3)
  real(dp) :: ux, uy, rx, ry
  real(dp) :: uqx(5,5), uqy(5,5), delta(2)
  real(dp) :: vk(2)
  real(dp), allocatable :: uq(:,:)

  Gn(1,:)=g1(:)
  Gn(2,:)=g12(:)
  Gn(3,:)=g12(:)-g1(:)
  Gn(4,:)=-Gn(1,:)
  Gn(5,:)=-Gn(2,:)
  Gn(6,:)=-Gn(3,:)    
  Gn(7,:)=Gn(1,:)
  Gn(8,:)=Gn(2,:)

  allocate(uq(ndim,2))
 
  Coordsnr(:,:) = Coords(:,:) 
 
  uqx=0.0_dp
  uqy=0.0_dp 

  ! Fourier components of the atomic displacement profile
  if(ntheta.EQ.2)then
  uqx(2,1)=-0.0007312154977677281
  uqy(2,1)=0.0004714203951774816
  uqx(3,1)=-1.154897989107238e-06
  uqy(3,1)=7.447553144970968e-07
  uqx(3,2)=-1.586938089071483e-06
  uqy(3,2)=0.0_dp
  uqx(4,1)=-2.434676944887857e-09
  uqy(4,1)=1.570103901518326e-09
  uqx(4,2)=-4.341715993276229e-09
  uqy(4,2)=1.000979455471334e-09
  uqx(4,3)=-4.341715993276229e-09
  uqy(4,3)=-1.000979455471334e-09
  uqx(5,1)=-5.701274058766884e-12
  uqy(5,1)=3.740789872160164e-12
  uqx(5,2)=-1.288484816335295e-11
  uqy(5,2)=4.352172282562038e-12
  uqx(5,3)=-1.480922507260647e-11
  uqy(5,3)=0.0_dp
  uqx(5,4)=-1.288484816335295e-11
  uqy(5,4)=-4.352172282562038e-12
  endif
  if(ntheta.EQ.3)then
  uqx(2,1)=-0.001435896542049262
  uqy(2,1)=0.0008971722346957707
  uqx(3,1)=-4.431771204029964e-06
  uqy(3,1)=2.770018809111695e-06
  uqx(3,2)=-5.980337522539451e-06
  uqy(3,2)=0.0_dp
  uqx(4,1)=-1.828082429017813e-08
  uqy(4,1)=1.142651350640817e-08
  uqx(4,2)=-3.192212955644209e-08
  uqy(4,2)=7.41532195221355e-09
  uqx(4,3)=-3.192212955644209e-08
  uqy(4,3)=-7.41532195221355e-09
  uqx(5,1)=-8.489859929513374e-11
  uqy(5,1)=5.318288726255601e-11
  uqx(5,2)=-1.864781523972493e-10
  uqy(5,2)=6.321663026418908e-11
  uqx(5,3)=-2.126102770742755e-10
  uqy(5,3)=0.0_dp
  uqx(5,4)=-1.864781523972493e-10
  uqy(5,4)=-6.321663026418908e-11
  endif
  if(ntheta.EQ.4)then
  uqx(2,1)=-0.002370539233839467
  uqy(2,1)=0.001455154743491223
  uqx(3,1)=-1.207488833535449e-05
  uqy(3,1)=7.415571455420628e-06
  uqx(3,2)=-1.606356439211949e-05
  uqy(3,2)=0.0_dp
  uqx(4,1)=-8.232092034607089e-08
  uqy(4,1)=5.05576290960723e-08
  uqx(4,2)=-1.412786099990225e-07
  uqy(4,2)=3.296261058695678e-08
  uqx(4,3)=-1.412786099990225e-07
  uqy(4,3)=-3.296261058695678e-08
  uqx(5,1)=-6.338434146746167e-10
  uqy(5,1)=3.894401693009466e-10
  uqx(5,2)=-1.363506478716756e-09
  uqy(5,2)=4.635827193166175e-10
  uqx(5,3)=-1.545001808355566e-09
  uqy(5,3)=0.0_dp
  uqx(5,4)=-1.363506478716756e-09
  uqy(5,4)=-4.635827193166175e-10
  endif
  if(ntheta.EQ.5)then
  uqx(2,1)=-0.00352744756912305
  uqy(2,1)=0.002140730764691146
  uqx(3,1)=-2.68010445377339e-05
  uqy(3,1)=1.627423402644026e-05
  uqx(3,2)=-3.516363993145865e-05
  uqy(3,2)=0.0_dp
  uqx(4,1)=-2.729706294022432e-07
  uqy(4,1)=1.657616285518053e-07
  uqx(4,2)=-4.604302998995516e-07
  uqy(4,2)=1.078115792758269e-07
  uqx(4,3)=-4.604302998995516e-07
  uqy(4,3)=-1.078115792758269e-07
  uqx(5,1)=-3.147029295689003e-09
  uqy(5,1)=1.911153328301163e-09
  uqx(5,2)=-6.629538953061466e-09
  uqy(5,2)=2.259053162518054e-09
  uqx(5,3)=-7.46551582954267e-09
  uqy(5,3)=0.0_dp
  uqx(5,4)=-6.629538953061466e-09
  uqy(5,4)=-2.259053162518054e-09
  endif
  if(ntheta.EQ.6)then
  uqx(2,1)=-0.004895854840040079
  uqy(2,1)=0.002947504859599212
  uqx(3,1)=-5.185833778140336e-05
  uqy(3,1)=3.124215464415819e-05
  uqx(3,2)=-6.704921137114243e-05
  uqy(3,2)=0.0_dp
  uqx(4,1)=-7.375225497266297e-07
  uqy(4,1)=4.443469426422899e-07
  uqx(4,2)=-1.22125720193429e-06
  uqy(4,2)=2.869558314022102e-07
  uqx(4,3)=-1.22125720193429e-06
  uqy(4,3)=-2.869558314022102e-07
  uqx(5,1)=-1.189937627779894e-08
  uqy(5,1)=7.168907495586338e-09
  uqx(5,2)=-2.451139179182853e-08
  uqy(5,2)=8.36909403937007e-09
  uqx(5,3)=-2.741742519309818e-08
  uqy(5,3)=0.0_dp
  uqx(5,4)=-2.451139179182853e-08
  uqy(5,4)=-8.36909403937007e-09
  endif
  if(ntheta.EQ.7)then
  uqx(2,1)=-0.006461449348026135
  uqy(2,1)=0.00386704590013867
  uqx(3,1)=-9.088953324281387e-05
  uqy(3,1)=5.443895214647056e-05
  uqx(3,2)=-0.000115650920643144
  uqy(3,2)=0.0_dp
  uqx(4,1)=-1.719155170137875e-06
  uqy(4,1)=1.029772845664776e-06
  uqx(4,2)=-2.790313685737744e-06
  uqy(4,2)=6.580010547012858e-07
  uqx(4,3)=-2.790313685737744e-06
  uqy(4,3)=-6.580010547012858e-07
  uqx(5,1)=-3.697520081119331e-08
  uqy(5,1)=2.214633072942774e-08
  uqx(5,2)=-7.434248153831823e-08
  uqy(5,2)=2.543295223277813e-08
  uqx(5,3)=-8.254619477080724e-08
  uqy(5,3)=0.0_dp
  uqx(5,4)=-7.434248153831823e-08
  uqy(5,4)=-2.543295223277813e-08
  endif
  if(ntheta.EQ.8)then
  uqx(2,1)=-0.008206072919252466
  uqy(2,1)=0.004888718057559328
  uqx(3,1)=-0.0001477464619472436
  uqy(3,1)=8.80996557824374e-05
  uqx(3,2)=-0.000184730257310846
  uqy(3,2)=0.0_dp
  uqx(4,1)=-3.581496930090418e-06
  uqy(4,1)=2.135790438804529e-06
  uqx(4,2)=-5.688326741335424e-06
  uqy(4,2)=1.34655485744528e-06
  uqx(4,3)=-5.688326741335424e-06
  uqy(4,3)=-1.34655485744528e-06
  uqx(5,1)=-9.895113648508195e-08
  uqy(5,1)=5.900265937736671e-08
  uqx(5,2)=-1.938529162557286e-07
  uqy(5,2)=6.645103335141035e-08
  uqx(5,3)=-2.135258580530001e-07
  uqy(5,3)=0.0_dp
  uqx(5,4)=-1.938529162557286e-07
  uqy(5,4)=-6.645103335141035e-08
  endif
  if(ntheta.EQ.9)then
  uqx(2,1)=-0.01010772196045918
  uqy(2,1)=0.005999689628423322
  uqx(3,1)=-0.0002262612975656868
  uqy(3,1)=0.0001344419648152459
  uqx(3,2)=-0.0002775244610317592
  uqy(3,2)=0.0_dp
  uqx(4,1)=-6.825276745678818e-06
  uqy(4,1)=4.055911695604714e-06
  uqx(4,2)=-1.059052135909834e-05
  uqy(4,2)=2.517294978569754e-06
  uqx(4,3)=-1.059052135909834e-05
  uqy(4,3)=-2.517294978569754e-06
  uqx(5,1)=-2.352002772937919e-07
  uqy(5,1)=1.397520149575992e-07
  uqx(5,2)=-4.48274468435581e-07
  uqy(5,2)=1.539852461003042e-07
  uqx(5,3)=-4.895322080001663e-07
  uqy(5,3)=0.0_dp
  uqx(5,4)=-4.48274468435581e-07
  uqy(5,4)=-1.539852461003042e-07
  endif
  if(ntheta.EQ.10)then
  uqx(2,1)=-0.01214094723296179
  uqy(2,1)=0.007185172779869293
  uqx(3,1)=-0.0003299929115003222
  uqy(3,1)=0.000195519493609438
  uqx(3,2)=-0.0003964101199496858
  uqy(3,2)=0.0_dp
  uqx(4,1)=-1.209151294434386e-05
  uqy(4,1)=7.165004163245931e-06
  uqx(4,2)=-1.830230180218659e-05
  uqy(4,2)=4.369241246891779e-06
  uqx(4,3)=-1.830230180218659e-05
  uqy(4,3)=-4.369241246891779e-06
  uqx(5,1)=-5.072472977662023e-07
  uqy(5,1)=3.005414175984017e-07
  uqx(5,2)=-9.393528199777028e-07
  uqy(5,2)=3.233772873927477e-07
  uqx(5,3)=-1.016494567114082e-06
  uqy(5,3)=0.0_dp
  uqx(5,4)=-9.393528199777028e-07
  uqy(5,4)=-3.233772873927477e-07
  endif
  if(ntheta.EQ.11)then
  uqx(2,1)=-0.01427768000921253
  uqy(2,1)=0.008428908168456773
  uqx(3,1)=-0.0004619778319639407
  uqy(3,1)=0.0002730781672413903
  uqx(3,2)=-0.0005426360143594519
  uqy(3,2)=0.0_dp
  uqx(4,1)=-2.014728185682702e-05
  uqy(4,1)=1.191074841386241e-05
  uqx(4,2)=-2.970961882633886e-05
  uqy(4,2)=7.124796626544694e-06
  uqx(4,3)=-2.970961882633886e-05
  uqy(4,3)=-7.124796626544694e-06
  uqx(5,1)=-1.008004944770031e-06
  uqy(5,1)=5.958424731241988e-07
  uqx(5,2)=-1.811980725955097e-06
  uqy(5,2)=6.251978084599946e-07
  uqx(5,3)=-1.942211173240885e-06
  uqy(5,3)=0.0_dp
  uqx(5,4)=-1.811980725955097e-06
  uqy(5,4)=-6.251978084599946e-07
  endif
  if(ntheta.EQ.12)then
  uqx(2,1)=-0.01648841976833508
  uqy(2,1)=0.00971385525211049
  uqx(3,1)=-0.0006245208306885359
  uqy(3,1)=0.0003684360208645578
  uqx(3,2)=-0.0007161707066312485
  uqy(3,2)=0.0_dp
  uqx(4,1)=-3.185339267975134e-05
  uqy(4,1)=1.879461360555224e-05
  uqx(4,2)=-4.570953565691623e-05
  uqy(4,2)=1.101356078295402e-05
  uqx(4,3)=-4.570953565691623e-05
  uqy(4,3)=-1.101356078295402e-05
  uqx(5,1)=-1.86730786016839e-06
  uqy(5,1)=1.101634412158447e-06
  uqx(5,2)=-3.256105072314334e-06
  uqy(5,2)=1.12609910984101e-06
  uqx(5,3)=-3.456108802116457e-06
  uqy(5,3)=0.0_dp
  uqx(5,4)=-3.256105072314334e-06
  uqy(5,4)=-1.12609910984101e-06
  endif
  if(ntheta.EQ.13)then
  uqx(2,1)=-0.01874363020388613
  uqy(2,1)=0.01102299870360727
  uqx(3,1)=-0.000819056300796783
  uqy(3,1)=0.0004824043177345613
  uqx(3,2)=-0.0009156910682583513
  uqy(3,2)=0.0_dp
  uqx(4,1)=-4.81174507641687e-05
  uqy(4,1)=2.834441813303126e-05
  uqx(4,2)=-6.713199189077586e-05
  uqy(4,2)=1.62532984542057e-05
  uqx(4,3)=-6.713199189077586e-05
  uqy(4,3)=-1.62532984542057e-05
  uqx(5,1)=-3.254182136977055e-06
  uqy(5,1)=1.916674289995359e-06
  uqx(5,2)=-5.502334381872796e-06
  uqy(5,2)=1.907490988568425e-06
  uqx(5,3)=-5.782489009202849e-06
  uqy(5,3)=0.0_dp
  uqx(5,4)=-5.502334381872796e-06
  uqy(5,4)=-1.907490988568425e-06
  endif
  if(ntheta.EQ.14)then
  uqx(2,1)=-0.02101513878443223
  uqy(2,1)=0.01234015250707122
  uqx(3,1)=-0.001046098478254697
  uqy(3,1)=0.0006152602359553973
  uqx(3,2)=-0.001138709526361856
  uqy(3,2)=0.0_dp
  uqx(4,1)=-6.983936154117285e-05
  uqy(4,1)=4.108260075599009e-05
  uqx(4,2)=-9.46663508367406e-05
  uqy(4,2)=2.303123144864916e-05
  uqx(4,3)=-9.46663508367406e-05
  uqy(4,3)=-2.303123144864916e-05
  uqx(5,1)=-5.374846134548899e-06
  uqy(5,1)=3.161278639182164e-06
  uqx(5,2)=-8.811118519862508e-06
  uqy(5,2)=3.06194712722794e-06
  uqx(5,3)=-9.167731209069543e-06
  uqy(5,3)=0.0_dp
  uqx(5,4)=-8.811118519862508e-06
  uqy(5,4)=-3.06194712722794e-06
  endif
  if(ntheta.EQ.15)then
  uqx(2,1)=-0.0232773402036116
  uqy(2,1)=0.01365064721466022
  uqx(3,1)=-0.001305280475347721
  uqy(3,1)=0.0007667708872687428
  uqx(3,2)=-0.001381810113794881
  uqy(3,2)=0.0_dp
  uqx(4,1)=-9.785794354653753e-05
  uqy(4,1)=5.749524297849592e-05
  uqx(4,2)=-0.0001288049100968271
  uqy(4,2)=3.148869813455972e-05
  uqx(4,3)=-0.0001288049100968271
  uqy(4,3)=-3.148869813455972e-05
  uqx(5,1)=-8.46637610884184e-06
  uqy(5,1)=4.973588375291042e-06
  uqx(5,2)=-1.345696760939882e-05
  uqy(5,2)=4.687804727529729e-06
  uqx(5,3)=-1.386349319225168e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-1.345696760939882e-05
  uqy(5,4)=-4.687804727529729e-06
  endif
  if(ntheta.EQ.16)then
  uqx(2,1)=-0.02550806137259484
  uqy(2,1)=0.01494181959537789
  uqx(3,1)=-0.001595465429771206
  uqy(3,1)=0.000936258825655531
  uqx(3,2)=-0.001640945873196756
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0001329064560728502
  uqy(4,1)=7.80063592484412e-05
  uqx(4,2)=-0.0001698108360147586
  uqy(4,2)=4.171126842036178e-05
  uqx(4,3)=-0.0001698108360147586
  uqy(4,3)=-4.171126842036178e-05
  uqx(5,1)=-1.278694321771364e-05
  uqy(5,1)=7.503853321298472e-06
  uqx(5,2)=-1.971030956794678e-05
  uqy(5,2)=6.882830706202828e-06
  uqx(5,3)=-2.010885042278717e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-1.971030956794678e-05
  uqy(5,4)=-6.882830706202828e-06
  endif
  if(ntheta.EQ.17)then
  uqx(2,1)=-0.02768902996623525
  uqy(2,1)=0.01620327280964095
  uqx(3,1)=-0.001914903012504294
  uqy(3,1)=0.001122693458825686
  uqx(3,2)=-0.001911747148401908
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0001755820633807187
  uqy(4,1)=0.0001029603397676965
  uqx(4,2)=-0.0002177119608478438
  uqy(4,2)=5.372503557832448e-05
  uqx(4,3)=-0.0002177119608478438
  uqy(4,3)=-5.372503557832448e-05
  uqx(5,1)=-1.860418596133933e-05
  uqy(5,1)=1.090765125005818e-05
  uqx(5,2)=-2.78198841681919e-05
  uqy(5,2)=9.737957708444016e-06
  uqx(5,3)=-2.811427798376512e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-2.78198841681919e-05
  uqy(5,4)=-9.737957708444016e-06
  endif
  if(ntheta.EQ.18)then
  uqx(2,1)=-0.0298059688900068
  uqy(2,1)=0.01742692140290434
  uqx(3,1)=-0.002261402797319784
  uqy(3,1)=0.001324791852622696
  uqx(3,2)=-0.00218979975575966
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0002263307798910572
  uqy(4,1)=0.0001326134069799172
  uqx(4,2)=-0.0002723167306780222
  uqy(4,2)=6.749853746454043e-05
  uqx(4,3)=-0.0002723167306780222
  uqy(4,3)=-6.749853746454043e-05
  uqx(5,1)=-2.618344775700557e-05
  uqy(5,1)=1.533904799393413e-05
  uqx(5,2)=-3.799806596259685e-05
  uqy(5,2)=1.333194246897647e-05
  uqx(5,3)=-3.804957109168239e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-3.799806596259685e-05
  uqy(5,4)=-1.333194246897647e-05
  endif
  if(ntheta.EQ.19)then
  uqx(2,1)=-0.03184839354005421
  uqy(2,1)=0.01860686639767005
  uqx(3,1)=-0.002632500871417905
  uqy(3,1)=0.001541115305242353
  uqx(3,2)=-0.00247086793286604
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0002854464623341651
  uqy(4,1)=0.0001671332222719948
  uqx(4,2)=-0.0003332456229776692
  uqy(4,2)=8.294893631494263e-05
  uqx(4,3)=-0.0003332456229776692
  uqy(4,3)=-8.294893631494263e-05
  uqx(5,1)=-3.5777331916109e-05
  uqy(5,1)=2.094453998211611e-05
  uqx(5,2)=-5.041053364039542e-05
  uqy(5,2)=1.772747608189362e-05
  uqx(5,3)=-5.003665347329104e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-5.041053364039542e-05
  uqy(5,4)=-1.772747608189362e-05
  endif
  if(ntheta.EQ.20)then
  uqx(2,1)=-0.03380921098354144
  uqy(2,1)=0.01973915809034029
  uqx(3,1)=-0.003025604178510975
  uqy(3,1)=0.001770152796829164
  uqx(3,2)=-0.002751053010400932
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0003530806004474008
  uqy(4,1)=0.0002066047359070089
  uqx(4,2)=-0.0003999707148415089
  uqy(4,2)=9.995079644234232e-05
  uqx(4,3)=-0.0003999707148415089
  uqy(4,3)=-9.995079644234232e-05
  uqx(5,1)=-4.761749968378625e-05
  uqy(5,1)=2.785831220463549e-05
  uqx(5,2)=-6.517066822235473e-05
  uqy(5,2)=2.296892009302885e-05
  uqx(5,3)=-6.414716777612842e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-6.517066822235473e-05
  uqy(5,4)=-2.296892009302885e-05
  endif
  if(ntheta.EQ.21)then
  uqx(2,1)=-0.03568421614472853
  uqy(2,1)=0.0208215014107083
  uqx(3,1)=-0.003438105278930305
  uqy(3,1)=0.002010387168015211
  uqx(3,2)=-0.003026890903307585
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0004292590682540321
  uqy(4,1)=0.0002510400416038367
  uqx(4,2)=-0.0004718571203498731
  uqy(4,2)=0.0001183459294007554
  uqx(4,3)=-0.0004718571203498731
  uqy(4,3)=-0.0001183459294007554
  uqx(5,1)=-6.190907493630538e-05
  uqy(5,1)=3.61990173774744e-05
  uqx(5,2)=-8.233825903866238e-05
  uqy(5,2)=2.908155778293243e-05
  uqx(5,3)=-8.040403158689358e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-8.233825903866238e-05
  uqy(5,4)=-2.908155778293243e-05
  endif
  if(ntheta.EQ.22)then
      uqx(2,1)= -3.7471561538033285E-002
      uqy(2,1)=  2.1852947132718681E-002
      uqx(3,1)= -3.8674671316928813E-003
      uqy(3,1)=  2.2603438961092790E-003
      uqx(3,2)= -3.2953989304928908E-003
      uqy(3,2)=  0.0000000000000000
      uqx(4,1)= -5.1390286043579824E-004
      uqy(4,1)=  3.0039050497730638E-004
      uqx(4,2)= -5.4820226802045503E-004
      uqy(4,2)=  1.3795319530098593E-004
      uqx(4,3)= -5.4820226802045503E-004
      uqy(4,3)= -1.3795319530098593E-004
      uqx(5,1)= -7.8827985776783088E-005
      uqy(5,1)=  4.6068260392008039E-005
      uqx(5,2)= -1.0192201661219018E-004
      uqy(5,2)=  3.6072174770764229E-005
      uqx(5,3)= -9.8786274203805616E-005
      uqy(5,3)=  0.0000000000000000
      uqx(5,4)= -1.0192201661219018E-004
      uqy(5,4)= -3.6072174770764229E-005
  endif
  if(ntheta.EQ.23)then
      uqx(2,1)= -3.9171248613853675E-002
      uqy(2,1)=  2.2833599959186702E-002
      uqx(3,1)= -4.3112776832871375E-003
      uqy(3,1)=  2.5186222740660215E-003
      uqx(3,2)= -3.5540820075879590E-003
      uqy(3,2)=  0.0000000000000000
      uqx(4,1)= -6.0684674564120337E-004
      uqy(4,1)=  3.5455765666248949E-004
      uqx(4,2)= -6.2826779626449060E-004
      uqy(4,2)=  1.5857712883648905E-004
      uqx(4,3)= -6.2826779626449060E-004
      uqy(4,3)= -1.5857712883648905E-004
      uqx(5,1)= -9.8518381717037002E-005
      uqy(5,1)=  5.7549139520913602E-005
      uqx(5,2)= -1.2388256346061724E-004
      uqy(5,2)=  4.3930047790534708E-005
      uqx(5,3)= -1.1923331918165264E-004
      uqy(5,3)=  0.0000000000000000
      uqx(5,4)= -1.2388256346061724E-004
      uqy(5,4)= -4.3930047790534708E-005
  endif
  if(ntheta.EQ.24)then
      uqx(2,1)= -4.0784680411922937E-002
      uqy(2,1)=  2.3764356805771322E-002
      uqx(3,1)= -4.7672873741273891E-003
      uqy(3,1)=  2.7839169205840812E-003
      uqx(3,2)= -3.8009168860933772E-003
      uqy(3,2)=  0.0000000000000000
      uqx(4,1)= -7.0786230269381919E-004
      uqy(4,1)=  4.1340658441093225E-004
      uqx(4,2)= -7.1130967374835976E-004
      uqy(4,2)=  1.8001594920585869E-004
      uqx(4,3)= -7.1130967374835976E-004
      uqy(4,3)= -1.8001594920585869E-004
      uqx(5,1)= -1.2109581826802820E-004
      uqy(5,1)=  7.0708105816613394E-005
      uqx(5,2)= -1.4814071815562693E-004
      uqy(5,2)=  5.2629469292820879E-005
      uqx(5,3)= -1.4165451131498082E-004
      uqy(5,3)=  0.0000000000000000
      uqx(5,4)= -1.4814071815562693E-004
      uqy(5,4)= -5.2629469292820879E-005
  endif
  if(ntheta.EQ.25)then
      uqx(2,1)= -4.2314275503676363E-002
      uqy(2,1)=  2.4646687206293511E-002
      uqx(3,1)= -5.2334239164346236E-003
      uqy(3,1)=  3.0550259837824490E-003
      uqx(3,2)= -4.0343157586948805E-003
      uqy(3,2)=  0.0000000000000000
      uqx(4,1)= -8.1667314921549940E-004
      uqy(4,1)=  4.7677479034404185E-004
      uqx(4,2)= -7.9659753583125406E-004
      uqy(4,2)=  2.0206753169306818E-004
      uqx(4,3)= -7.9659753583125406E-004
      uqy(4,3)= -2.0206753169306818E-004
      uqx(5,1)= -1.4664751584463312E-004
      uqy(5,1)=  8.5595146227575536E-005
      uqx(5,2)= -1.7458254500867196E-004
      uqy(5,2)=  6.2131570401939717E-005
      uqx(5,3)= -1.6593414597602474E-004
      uqy(5,3)=  0.0000000000000000
      uqx(5,4)= -1.7458254500867196E-004
      uqy(5,4)= -6.2131570401939717E-005
  endif
  if(ntheta.EQ.26)then
      uqx(2,1)= -4.3763154783815199E-002
      uqy(2,1)=  2.5482453046564175E-002
      uqx(3,1)= -5.7077974506492624E-003
      uqy(3,1)=  3.3308539124622160E-003
      uqx(3,2)= -4.2530836670724090E-003
      uqy(3,2)=  0.0000000000000000
      uqx(4,1)= -9.3297026737261708E-004
      uqy(4,1)=  5.4448108076819474E-004
      uqx(4,2)= -8.8343066005766308E-004
      uqy(4,2)=  2.2453437830041617E-004
      uqx(4,3)= -8.8343066005766308E-004
      uqy(4,3)= -2.2453437830041617E-004
      uqx(5,1)= -1.7523533791535451E-004
      uqy(5,1)=  1.0224552756457173E-004
      uqx(5,2)= -2.0306641416121418E-004
      uqy(5,2)=  7.2386717005680300E-005
      uqx(5,3)= -1.9193831138938561E-004
      uqy(5,3)=  0.0000000000000000
      uqx(5,4)= -2.0306641416121418E-004
      uqy(5,4)= -7.2386717005680300E-005
  endif
  if(ntheta.EQ.27)then
      uqx(2,1)= -4.5134892181384466E-002
      uqy(2,1)=  2.6273765959396265E-002
      uqx(3,1)= -6.1886967887041519E-003
      uqy(3,1)=  3.6104090982189683E-003
      uqx(3,2)= -4.4563723758183297E-003
      uqy(3,2)=  0.0000000000000000
      uqx(4,1)= -1.0564240178767453E-003
      uqy(4,1)=  6.1633252655025406E-004
      uqx(4,2)= -9.7114868671522236E-004
      uqy(4,2)=  2.4722732168313730E-004
      uqx(4,3)= -9.7114868671522236E-004
      uqy(4,3)= -2.4722732168313730E-004
      uqx(5,1)= -2.0689871665405404E-004
      uqy(5,1)=  1.2068150811888939E-004
      uqx(5,2)= -2.3342935306930540E-004
      uqy(5,2)=  8.3336740327999767E-005
      uqx(5,3)= -2.1952056758095347E-004
      uqy(5,3)=  0.0000000000000000
      uqx(5,4)= -2.3342935306930540E-004
      uqy(5,4)= -8.3336740327999767E-005
  endif
  if(ntheta.EQ.28)then
      uqx(2,1)= -4.6433321768224768E-002
      uqy(2,1)=  2.7022877302892028E-002
      uqx(3,1)= -6.6745804119781290E-003
      uqy(3,1)=  3.8927985245745503E-003
      uqx(3,2)= -4.6436346742308558E-003
      uqy(3,2)=  0.0000000000000000
      uqx(4,1)= -1.1866935904011658E-003
      uqy(4,1)=  6.9212993050260974E-004
      uqx(4,2)= -1.0591384787539894E-003
      uqy(4,2)=  2.6996820723435564E-004
      uqx(4,3)= -1.0591384787539894E-003
      uqy(4,3)= -2.6996820723435564E-004
      uqx(5,1)= -2.4165765845836861E-004
      uqy(5,1)=  1.4091409176613024E-004
      uqx(5,2)= -2.6549284977472415E-004
      uqy(5,2)=  9.4917013558184400E-005
      uqx(5,3)= -2.4852679837027429E-004
      uqy(5,3)=  0.0000000000000000
      uqx(5,4)= -2.6549284977472415E-004
      uqy(5,4)= -9.4917013558184400E-005
  endif
  if(ntheta.EQ.29)then
      uqx(2,1)= -4.7662392078997758E-002
      uqy(2,1)=  2.7732095343409147E-002
      uqx(3,1)= -7.1640645118234338E-003
      uqy(3,1)=  4.1772207415418164E-003
      uqx(3,2)= -4.8145811869661570E-003
      uqy(3,2)=  0.0000000000000000
      uqx(4,1)= -1.3234342016656580E-003
      uqy(4,1)=  7.7167198538604460E-004
      uqx(4,2)= -1.1468379861261022E-003
      uqy(4,2)=  2.9259173325616312E-004
      uqx(4,3)= -1.1468379861261022E-003
      uqy(4,3)= -2.9259173325616312E-004
      uqx(5,1)= -2.7951564618481208E-004
      uqy(5,1)=  1.6294471892770439E-004
      uqx(5,2)= -2.9906799566926711E-004
      uqy(5,2)=  1.0705832709851686E-004
      uqx(5,3)= -2.7879922322344216E-004
      uqy(5,3)=  0.0000000000000000
      uqx(5,4)= -2.9906799566926711E-004
      uqy(5,4)= -1.0705832709851686E-004
  endif
  if(ntheta.EQ.30)then
      uqx(2,1)= -4.8826058753127836E-002
      uqy(2,1)=  2.8403724470307239E-002
      uqx(3,1)= -7.6559096797667245E-003
      uqy(3,1)=  4.4629580893927974E-003
      uqx(3,2)= -4.9691407814230594E-003
      uqy(3,2)=  0.0000000000000000
      uqx(4,1)= -1.4663023959723578E-003
      uqy(4,1)=  8.5475833043351123E-004
      uqx(4,2)= -1.2337378598176337E-003
      uqy(4,2)=  3.1494661456099707E-004
      uqx(4,3)= -1.2337378598176337E-003
      uqy(4,3)= -3.1494661456099707E-004
      uqx(5,1)= -3.2046232246823515E-004
      uqy(5,1)=  1.8676682796613889E-004
      uqx(5,2)= -3.3395993507357336E-004
      uqy(5,2)=  1.1968854347660691E-004
      uqx(5,3)= -3.1017961799516920E-004
      uqy(5,3)=  0.0000000000000000
      uqx(5,4)= -3.3395993507357336E-004
      uqy(5,4)= -1.1968854347660691E-004
  endif
  if(ntheta.EQ.31)then
      uqx(2,1)= -4.9928207554533525E-002
      uqy(2,1)=  2.9040021838646290E-002
      uqx(3,1)= -8.1490073273388936E-003
      uqy(3,1)=  4.7493687904768891E-003
      uqx(3,2)= -5.1074250061246308E-003
      uqy(3,2)=  0.0000000000000000
      uqx(4,1)= -1.6149597985871669E-003
      uqy(4,1)=  9.4119171099780324E-004
      uqx(4,2)= -1.3193814253729618E-003
      uqy(4,2)=  3.3689621236779977E-004
      uqx(4,3)= -1.3193814253729618E-003
      uqy(4,3)= -3.3689621236779977E-004
      uqx(5,1)= -3.6447589010750524E-004
      uqy(5,1)=  2.1236725049130582E-004
      uqx(5,2)= -3.6997163938217169E-004
      uqy(5,2)=  1.3273403147615894E-004
      uqx(5,3)= -3.4251182702598322E-004
      uqy(5,3)=  0.0000000000000000
      uqx(5,4)= -3.6997163938217169E-004
      uqy(5,4)= -1.3273403147615894E-004
  endif
  if(ntheta.EQ.32)then
      uqx(2,1)= -5.0972601014766056E-002
      uqy(2,1)=  2.9643167537421276E-002
      uqx(3,1)= -8.6423665252589980E-003
      uqy(3,1)=  5.0358793026496528E-003
      uqx(3,2)= -5.2296965906485542E-003
      uqy(3,2)=  0.0000000000000000
      uqx(4,1)= -1.7690756422560100E-003
      uqy(4,1)=  1.0307794278985045E-003
      uqx(4,2)= -1.4033634968132719E-003
      uqy(4,2)=  3.5831874920469177E-004
      uqx(4,3)= -1.4033634968132719E-003
      uqy(4,3)= -3.5831874920469177E-004
      uqx(5,1)= -4.1152520287070926E-004
      uqy(5,1)=  2.3972742565373705E-004
      uqx(5,2)= -4.0690705252185934E-004
      uqy(5,2)=  1.4612089010502584E-004
      uqx(5,3)= -3.7564366300331148E-004
      uqy(5,3)=  0.0000000000000000
      uqx(5,4)= -4.0690705252185934E-004
      uqy(5,4)= -1.4612089010502584E-004
  endif
  if(ntheta.EQ.33)then
      uqx(2,1)= -5.1962843147254065E-002
      uqy(2,1)=  3.0215245083260955E-002
      uqx(3,1)= -9.1351016727408973E-003
      uqy(3,1)=  5.3219771671499990E-003
      uqx(3,2)= -5.3363418074061459E-003
      uqy(3,2)=  0.0000000000000000
      uqx(4,1)= -1.9283283457040433E-003
      uqy(4,1)=  1.1233342386024689E-003
      uqx(4,2)= -1.4853283997341403E-003
      uqy(4,2)=  3.7910720499739672E-004
      uqx(4,3)= -1.4853283997341403E-003
      uqy(4,3)= -3.7910720499739672E-004
      uqx(5,1)= -4.6157154577173112E-004
      uqy(5,1)=  2.6882443331406765E-004
      uqx(5,2)= -4.4457366870887026E-004
      uqy(5,2)=  1.5977597888436947E-004
      uqx(5,3)= -4.0942829263365720E-004
      uqy(5,3)=  0.0000000000000000
      uqx(5,4)= -4.4457366870887026E-004
      uqy(5,4)= -1.5977597888436947E-004
  endif
  if(ntheta.EQ.34)then
      uqx(2,1)=    -5.2902357777429854E-002
      uqy(2,1)=     3.0758229677107351E-002
      uqx(3,1)=    -9.6264212142102550E-003
      uqy(3,1)=     5.6072044725382752E-003
      uqx(3,2)=    -5.4278463797914674E-003
      uqy(3,2)=     0.0000000000000000
      uqx(4,1)=    -2.0924063788045214E-003
      uqy(4,1)=     1.2186748464810043E-003
      uqx(4,2)=    -1.5649674801196977E-003
      uqy(4,2)=     3.9916897105356127E-004
      uqx(4,3)=    -1.5649674801196977E-003
      uqy(4,3)=    -3.9916897105356127E-004
      uqx(5,1)=    -5.1457012054211047E-004
      uqy(5,1)=     2.9963185555462353E-004
      uqx(5,2)=    -4.8278460816578568E-004
      uqy(5,2)=     1.7362777338279618E-004
      uqx(5,3)=    -4.4372520064773737E-004
      uqy(5,3)=     0.0000000000000000
      uqx(5,4)=    -4.8278460816578568E-004
      uqy(5,4)=    -1.7362777338279618E-004
  endif
  if(ntheta.EQ.35)then
      uqx(2,1)=-5.3794376976085793E-002
      uqy(2,1)=3.1273982206297665E-002
      uqx(3,1)=-1.0115617485696799E-002
      uqy(3,1)=5.8911519789315529E-003
      uqx(3,2)=-5.5047745794039836E-003
      uqy(3,2)=0.0000000000000000
      uqx(4,1)=-2.2610086021982588E-003
      uqy(4,1)=1.3166260869920154E-003
      uqx(4,2)=-1.6420163024753152E-003
      uqy(4,2)=4.1842532338285868E-004
      uqx(4,3)=-1.6420163024753152E-003
      uqy(4,3)=-4.1842532338285868E-004
      uqx(5,1)=-5.7047125819495415E-004
      uqy(5,1)=3.3212047948659196E-004
      uqx(5,2)=-5.2136025511405142E-004
      uqy(5,2)=1.8760706555815432E-004
      uqx(5,3)=-4.7840081370311934E-004
      uqy(5,3)=0.0000000000000000
      uqx(5,4)=-5.2136025511405142E-004
      uqy(5,4)=-1.8760706555815432E-004
  endif
  if(ntheta.EQ.36)then
      uqx(2,1)=-5.4641936878380391E-002
      uqy(2,1)=3.1764247439854169E-002
      uqx(3,1)=-1.0602057725490647E-002
      uqy(3,1)=6.1734539194323035E-003
      uqx(3,2)=-5.5677511368476410E-003
      uqy(3,2)=0.0000000000000000
      uqx(4,1)=-2.4338442693324466E-003
      uqy(4,1)=1.4170189197121887E-003
      uqx(4,2)=-1.7162516835167235E-003
      uqy(4,2)=4.3681075786328644E-004
      uqx(4,3)=-1.7162516835167235E-003
      uqy(4,3)=-4.3681075786328644E-004
      uqx(5,1)=-6.2922141624780837E-004
      uqy(5,1)=3.6625887498432409E-004
      uqx(5,2)=-5.6012951988682709E-004
      uqy(5,2)=2.0164752506796602E-004
      uqx(5,3)=-5.1332886882180873E-004
      uqy(5,3)=0.0000000000000000
      uqx(5,4)=-5.6012951988682709E-004
      uqy(5,4)=-2.0164752506796602E-004
  endif
  if(ntheta.EQ.37)then
    uqx(2,1)=-5.5447878759301759E-002
    uqy(2,1)=3.2230655181018031E-002
    uqx(3,1)=-1.1085176099977573E-002
    uqy(3,1)=6.4537833913286312E-003
    uqx(3,2)=-5.6174456711525758E-003
    uqy(3,2)=0.0000000000000000
    uqx(4,1)=-2.6106327129279995E-003
    uqy(4,1)=1.5196902385646830E-003
    uqx(4,2)=-1.7874886689443466E-003
    uqy(4,2)=4.5427224122256848E-004
    uqx(4,3)=-1.7874886689443466E-003
    uqy(4,3)=-4.5427224122256848E-004
    uqx(5,1)=-6.9076391337295154E-004
    uqy(5,1)=4.0201381991562617E-004
    uqx(5,2)=-5.9893077476703364E-004
    uqy(5,2)=2.1568614505137565E-004
    uqx(5,3)=-5.4839055176817146E-004
    uqy(5,3)=0.0000000000000000
    uqx(5,4)=-5.9893077476703364E-004
    uqy(5,4)=-2.1568614505137565E-004
  endif
  if(ntheta.EQ.38)then
      uqx(2,1)=-5.6214853825552059E-002
      uqy(2,1)=3.2674723524759489E-002
      uqx(3,1)=-1.1564466825752056E-002
      uqy(3,1)=6.7318483805712615E-003
      uqx(3,2)=-5.6545592578638908E-003
      uqy(3,2)=0.0000000000000000
      uqx(4,1)=-2.7911029766045301E-003
      uqy(4,1)=1.6244826515550953E-003
      uqx(4,2)=-1.8555775218777036E-003
      uqy(4,2)=4.7076838115865015E-004
      uqx(4,3)=-1.8555775218777036E-003
      uqy(4,3)=-4.7076838115865015E-004
      uqx(5,1)=-7.5503956722783544E-004
      uqy(5,1)=4.3935066935865469E-004
      uqx(5,2)=-6.3761251978989718E-004
      uqy(5,2)=2.2966357654393945E-004
      uqx(5,3)=-5.8347451708535094E-004
      uqy(5,3)=0.0000000000000000
      uqx(5,4)=-6.3761251978989718E-004
      uqy(5,4)=-2.2966357654393945E-004
  endif
  if(ntheta.EQ.39)then
      uqx(2,1)=-5.6945330466828245E-002
      uqy(2,1)=3.3097863473680995E-002
      uqx(3,1)=-1.2039478126645067E-002
      uqy(3,1)=7.0073882687826375E-003
      uqx(3,2)=-5.6798129387907650E-003
      uqy(3,2)=0.0000000000000000
      uqx(4,1)=-2.9749932964843501E-003
      uqy(4,1)=1.7312441742406643E-003
      uqx(4,2)=-1.9204007778897718E-003
      uqy(4,2)=4.8626856663909642E-004
      uqx(4,3)=-1.9204007778897718E-003
      uqy(4,3)=-4.8626856663909642E-004
      uqx(5,1)=-8.2198713275973185E-004
      uqy(5,1)=4.7823360900544292E-004
      uqx(5,2)=-6.7603381259981592E-004
      uqy(5,2)=2.4352437605165144E-004
      uqx(5,3)=-6.1847677131471540E-004
      uqy(5,3)=0.0000000000000000
      uqx(5,4)=-6.7603381259981592E-004
      uqy(5,4)=-2.4352437605165144E-004
  endif
  if(ntheta.EQ.40)then
      uqx(2,1)=-5.7641603101368837E-002
      uqy(2,1)=3.3501384435461808E-002
      uqx(3,1)=-1.2509807019621438E-002
      uqy(3,1)=7.2801708174011770E-003
      uqx(3,2)=-5.6939378938862511E-003
      uqy(3,2)=0.0000000000000000
      uqx(4,1)=-3.1620505852235823E-003
      uqy(4,1)=1.8398279253301218E-003
      uqx(4,2)=-1.9818703975037046E-003
      uqy(4,2)=5.0075208191544451E-004
      uqx(4,3)=-1.9818703975037046E-003
      uqy(4,3)=-5.0075208191544451E-004
      uqx(5,1)=-8.9154364897433072E-004
      uqy(5,1)=5.1862585556785104E-004
      uqx(5,2)=-7.1406450376749575E-004
      uqy(5,2)=2.5721717160551692E-004
      uqx(5,3)=-6.5330048929207170E-004
      uqy(5,3)=0.0000000000000000
      uqx(5,4)=-7.1406450376749575E-004
      uqy(5,4)=-2.5721717160551692E-004
  endif

  !Update coordinates
  do i=1,ndim/2

    rx=Coords(i,1)
    ry=Coords(i,2)

    delta(1)=(0.d0,0.d0)
    delta(2)=(0.d0,0.d0)

    do h=1,5,2
      do n=2,5
        do m=1,n-1

          vk(1)=(n-1)*Gn(h,1)+(m-1)*Gn(h+2,1)
          vk(2)=(n-1)*Gn(h,2)+(m-1)*Gn(h+2,2)

          ux=cos((h-1)*pi/3.0_dp)*uqx(n,m)-sin((h-1)*pi/3.0_dp)*uqy(n,m)
          uy=sin((h-1)*pi/3.0_dp)*uqx(n,m)+cos((h-1)*pi/3.0_dp)*uqy(n,m)

          delta(1)=delta(1)-2*ux*sin(vk(1)*rx+vk(2)*ry)
          delta(2)=delta(2)-2*uy*sin(vk(1)*rx+vk(2)*ry)

        end do
      end do
    end do

    Coords(i,1)=rx+delta(1)/2
    Coords(i,2)=ry+delta(2)/2

    Coords(i+ndim/2,1)=-Coords(i,1)!rx-delta(1)/2
    Coords(i+ndim/2,2)=Coords(i,2) !ry-delta(2)/2
   
  enddo

  Nk=0
  do i=ndim/2+1,3*ndim/4
    do j=ndim/2+1,3*ndim/4
      if(sqrt((Coords(i,1) - Coordsnr(j,1))**2 + (Coords(i,2) - Coordsnr(j,2))**2).lt.0.25_dp)then
        Nk = Nk+1
        Coordstemp(j,:) = Coords(i,:)
        exit
      endif 
    enddo
  enddo
  
  Nk=Nk+1
  Coordstemp(3*ndim/4,1) = Coords(3*ndim/4,1)
  Coordstemp(3*ndim/4,2) = -Coords(3*ndim/4,2)
  Coordstemp(3*ndim/4,3) = Coords(3*ndim/4,3)
  
  if(Nk.ne.ndim/4)then
    write(*,*) 'ERROR RELAXATION', Nk
  endif

  Nk=0
  do i=3*ndim/4+1,ndim
    do j=3*ndim/4+1,ndim
      if(sqrt((Coords(i,1) - Coordsnr(j,1))**2 + (Coords(i,2) - Coordsnr(j,2))**2).lt.0.25_dp)then
        Nk = Nk+1
        Coordstemp(j,:) = Coords(i,:)
        exit
      endif 
    enddo
  enddo
  
  Nk=Nk+1
  Coordstemp(ndim,1) = Coords(ndim,1)
  Coordstemp(ndim,2) = -Coords(ndim,2)
  Coordstemp(ndim,3) = Coords(ndim,3)
  
  if(Nk.ne.ndim/4)then
    write(*,*) 'ERROR RELAXATION', Nk
  endif

  Coords(ndim/2+1:ndim,:) = Coordstemp(ndim/2+1:ndim,:)


end subroutine LatticeRelaxationKoshino  

subroutine LatticeRelaxationCarr(Coords,g1,g12)

  real(dp), intent(inout) :: Coords(ndim,3)
  real(dp), intent(in)    :: g1(2), g12(2) 

  real(dp) :: Gn(8,2)
  integer(dp) :: i,j, m, n, Nk, h
  real(dp) :: Coordsnr(ndim,3), Coordstemp(ndim,3)
  real(dp) :: ux, uy, rx, ry
  real(dp) :: uqx(5,5), uqy(5,5), delta(2)
  real(dp) :: vk(2)
  real(dp), allocatable :: uq(:,:)

  Gn(1,:)=g1(:)
  Gn(2,:)=g12(:)
  Gn(3,:)=g12(:)-g1(:)
  Gn(4,:)=-Gn(1,:)
  Gn(5,:)=-Gn(2,:)
  Gn(6,:)=-Gn(3,:)    
  Gn(7,:)=Gn(1,:)
  Gn(8,:)=Gn(2,:)

  allocate(uq(ndim,2))
 
  Coordsnr(:,:) = Coords(:,:) 
 
  uqx=0.0_dp
  uqy=0.0_dp 

  ! Fourier components of the atomic displacement profile
 if (ntheta.eq.          20 )then
 uqx(2,1)=  -1.2330860177910182E-002
 uqy(2,1)=   7.2063554058462605E-003
 uqx(3,1)=  -2.0660797221745415E-004
 uqy(3,1)=   1.2137941856539299E-004
 uqx(3,2)=   3.0817470087879243E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=   2.3292649044767669E-006
 uqy(4,1)=  -1.3297080603197937E-006
 uqx(4,2)=   4.5899567069562079E-005
 uqy(4,2)=  -6.4738474373867059E-006
 uqx(4,3)=   4.5899567069562079E-005
 uqy(4,3)=   6.4738474373867059E-006
 uqx(5,1)=   5.3181324805917406E-007
 uqy(5,1)=  -3.1012009191337748E-007
 uqx(5,2)=   3.2347840817742606E-006
 uqy(5,2)=  -7.6501142398045918E-007
 uqx(5,3)=   4.6755639825218119E-006
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   3.2347840817742606E-006
 uqy(5,4)=   7.6501142398045918E-007
 endif
 if (ntheta.eq.          21 )then
 uqx(2,1)=  -1.3477849743181924E-002
 uqy(2,1)=   7.8716219360763628E-003
 uqx(3,1)=  -2.6534914964265814E-004
 uqy(3,1)=   1.5571038787831639E-004
 uqx(3,2)=   2.9290950239227112E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=   1.2437842505030771E-006
 uqy(4,1)=  -6.8514794574154778E-007
 uqx(4,2)=   5.2291319607017092E-005
 uqy(4,2)=  -7.1548245433390379E-006
 uqx(4,3)=   5.2291319607017092E-005
 uqy(4,3)=   7.1548245433390379E-006
 uqx(5,1)=   6.5946594105574240E-007
 uqy(5,1)=  -3.8393449878868099E-007
 uqx(5,2)=   4.1782449759769181E-006
 uqy(5,2)=  -9.7971455658410332E-007
 uqx(5,3)=   6.1403924752653818E-006
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   4.1782449759769181E-006
 uqy(5,4)=   9.7971455658410332E-007
 endif
 if (ntheta.eq.          22 )then
 uqx(2,1)=  -1.4656130951246866E-002
 uqy(2,1)=   8.5547786937535050E-003
 uqx(3,1)=  -3.3359290900504313E-004
 uqy(3,1)=   1.9556471781458039E-004
 uqx(3,2)=   2.6977852641566136E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -4.8212285048685738E-007
 uqy(4,1)=   3.3439106098250644E-007
 uqx(4,2)=   5.8723818108983004E-005
 uqy(4,2)=  -7.7512104506240507E-006
 uqx(4,3)=   5.8723818108983004E-005
 uqy(4,3)=   7.7512104506240507E-006
 uqx(5,1)=   7.9315017935500992E-007
 uqy(5,1)=  -4.6097877788023966E-007
 uqx(5,2)=   5.2908212876999923E-006
 uqy(5,2)=  -1.2276873601027506E-006
 uqx(5,3)=   7.8952590917890556E-006
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   5.2908212876999923E-006
 uqy(5,4)=   1.2276873601027506E-006
 endif
 if (ntheta.eq.          23 )then
 uqx(2,1)=  -1.5861326738265374E-002
 uqy(2,1)=   9.2532941694343316E-003
 uqx(3,1)=  -4.1194485132302084E-004
 uqy(3,1)=   2.4129133975456986E-004
 uqx(3,2)=   2.3840253137431738E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -3.0106327238015306E-006
 uqy(4,1)=   1.8232330144996248E-006
 uqx(4,2)=   6.5027319824044743E-005
 uqy(4,2)=  -8.2216916685826734E-006
 uqx(4,3)=   6.5027319824044743E-005
 uqy(4,3)=   8.2216916685826734E-006
 uqx(5,1)=   9.2298837597126649E-007
 uqy(5,1)=  -5.3543973083821619E-007
 uqx(5,2)=   6.5766129828380905E-006
 uqy(5,2)=  -1.5071360675399520E-006
 uqx(5,3)=   9.9577245321376886E-006
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   6.5766129828380905E-006
 uqy(5,4)=   1.5071360675399520E-006
 endif
 if (ntheta.eq.          24 )then
 uqx(2,1)=  -1.7088878807244690E-002
 uqy(2,1)=   9.9645338985212360E-003
 uqx(3,1)=  -5.0093710421777982E-004
 uqy(3,1)=   2.9319637852414482E-004
 uqx(3,2)=   1.9852752891428142E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -6.5202395608067339E-006
 uqy(4,1)=   3.8851039738690907E-006
 uqx(4,2)=   7.1013883635102147E-005
 uqy(4,2)=  -8.5218958519750018E-006
 uqx(4,3)=   7.1013883635102147E-005
 uqy(4,3)=   8.5218958519750018E-006
 uqx(5,1)=   1.0347485144215188E-006
 uqy(5,1)=  -5.9897051517167228E-007
 uqx(5,2)=   8.0334446288719718E-006
 uqy(5,2)=  -1.8141951388210562E-006
 uqx(5,3)=   1.2338033006189197E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   8.0334446288719718E-006
 uqy(5,4)=   1.8141951388210562E-006
 endif
 if (ntheta.eq.          25 )then
 uqx(2,1)=  -1.8334117835630912E-002
 uqy(2,1)=   1.0685801112939312E-002
 uqx(3,1)=  -6.0101602548890310E-004
 uqy(3,1)=   3.5153607162012162E-004
 uqx(3,2)=   1.5003251512483640E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -1.1203313714493506E-005
 uqy(4,1)=   6.6317042211664737E-006
 uqx(4,2)=   7.6482076940620112E-005
 uqy(4,2)=  -8.6055436668371758E-006
 uqx(4,3)=   7.6482076940620112E-005
 uqy(4,3)=   8.6055436668371758E-006
 uqx(5,1)=   1.1091826249022040E-006
 uqy(5,1)=  -6.4030967312794144E-007
 uqx(5,2)=   9.6519161972274835E-006
 uqy(5,2)=  -2.1426999049442087E-006
 uqx(5,3)=   1.5037773449454520E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   9.6519161972274835E-006
 uqy(5,4)=   2.1426999049442087E-006
 endif
 if (ntheta.eq.          26 )then
 uqx(2,1)=  -1.9592339174721231E-002
 uqy(2,1)=   1.1414380505389916E-002
 uqx(3,1)=  -7.1253240889183227E-004
 uqy(3,1)=   4.1651114507332022E-004
 uqx(3,2)=   9.2932651708509557E-005
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -1.7262965830187229E-005
 uqy(4,1)=   1.0180880421053244E-005
 uqx(4,2)=   8.1222451029173295E-005
 uqy(4,2)=  -8.4257436652099532E-006
 uqx(4,3)=   8.1222451029173295E-005
 uqy(4,3)=   8.4257436652099532E-006
 uqx(5,1)=   1.1214801345017206E-006
 uqy(5,1)=  -6.4496745793858162E-007
 uqx(5,2)=   1.1414733275946095E-005
 uqy(5,2)=  -2.4840537198367436E-006
 uqx(5,3)=   1.8048878075634542E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   1.1414733275946095E-005
 uqy(5,4)=   2.4840537198367436E-006
 endif
 if (ntheta.eq.          27 )then
 uqx(2,1)=  -2.0858880297683104E-002
 uqy(2,1)=   1.2147582940679243E-002
 uqx(3,1)=  -8.3573470555703465E-004
 uqy(3,1)=   4.8826294096080655E-004
 uqx(3,2)=   2.7377702450722934E-005
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -2.4909356928400127E-005
 uqy(4,1)=   1.4654479613094008E-005
 uqx(4,2)=   8.5023481689503201E-005
 uqy(4,2)=  -7.9363587287342227E-006
 uqx(4,3)=   8.5023481689503201E-005
 uqy(4,3)=   7.9363587287342227E-006
 uqx(5,1)=   1.0408844747303074E-006
 uqy(5,1)=  -5.9500734155654360E-007
 uqx(5,2)=   1.3296381750714867E-005
 uqy(5,2)=  -2.8272060961681220E-006
 uqx(5,3)=   2.1353040939838830E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   1.3296381750714867E-005
 uqy(5,4)=   2.8272060961681220E-006
 endif
 if (ntheta.eq.          28 )then
 uqx(2,1)=  -2.2129196139014162E-002
 uqy(2,1)=   1.2882788887369498E-002
 uqx(3,1)=  -9.7076556394750475E-004
 uqy(3,1)=   5.6687146925068382E-004
 uqx(3,2)=  -4.6354130209459111E-005
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -3.4355640549775946E-005
 uqy(4,1)=   2.0175993364914671E-005
 uqx(4,2)=   8.7677634149330643E-005
 uqy(4,2)=  -7.0933662552843666E-006
 uqx(4,3)=   8.7677634149330643E-005
 uqy(4,3)=   7.0933662552843666E-006
 uqx(5,1)=   8.3051131376083082E-007
 uqy(5,1)=  -4.6894472311642965E-007
 uqx(5,2)=   1.5263184777237425E-005
 uqy(5,2)=  -3.1587494164919436E-006
 uqx(5,3)=   2.4921605001521162E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   1.5263184777237425E-005
 uqy(5,4)=   3.1587494164919436E-006
 endif
 if (ntheta.eq.          29 )then
 uqx(2,1)=  -2.3398928719575527E-002
 uqy(2,1)=   1.3617488492555174E-002
 uqx(3,1)=  -1.1176617546509483E-003
 uqy(3,1)=   6.5235541951759479E-004
 uqx(3,2)=  -1.2786650353677659E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -4.5813750763561952E-005
 uqy(4,1)=   2.6868116195381131E-005
 uqx(4,2)=   8.8987208387438240E-005
 uqy(4,2)=  -5.8561357282975438E-006
 uqx(4,3)=   8.8987208387438240E-005
 uqy(4,3)=   5.8561357282975438E-006
 uqx(5,1)=   4.4739209562586133E-007
 uqy(5,1)=  -2.4177630899923655E-007
 uqx(5,2)=   1.7273748095154128E-005
 uqy(5,2)=  -3.4631319300386581E-006
 uqx(5,3)=   2.8715926765090612E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   1.7273748095154128E-005
 uqy(5,4)=   3.4631319300386581E-006
 endif
 if (ntheta.eq.          30 )then
 uqx(2,1)=  -2.4663968024620533E-002
 uqy(2,1)=   1.4349316557387061E-002
 uqx(3,1)=  -1.2763573132133310E-003
 uqy(3,1)=   7.4467403323366990E-004
 uqx(3,2)=  -2.1665864038085869E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -5.9490255992798609E-005
 uqy(4,1)=   3.4850345661711466E-005
 uqx(4,2)=   8.8769651796954534E-005
 uqy(4,2)=  -4.1885562358857108E-006
 uqx(4,3)=   8.8769651796954534E-005
 uqy(4,3)=   4.1885562358857108E-006
 uqx(5,1)=  -1.5725087317644755E-007
 uqy(5,1)=   1.1485648259070322E-007
 uqx(5,2)=   1.9279767732597909E-005
 uqy(5,2)=  -3.7229752040447505E-006
 uqx(5,3)=   3.2688188304589006E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   1.9279767732597909E-005
 uqy(5,4)=   3.7229752040447505E-006
 endif
 if (ntheta.eq.          31 )then
 uqx(2,1)=  -2.5920501908694788E-002
 uqy(2,1)=   1.5076081139368567E-002
 uqx(3,1)=  -1.4466895320886283E-003
 uqy(3,1)=   8.4373062067896558E-004
 uqx(3,2)=  -3.1214153806844840E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -7.5582483455627225E-005
 uqy(4,1)=   4.4236742464516627E-005
 uqx(4,2)=   8.6862085239109360E-005
 uqy(4,2)=  -2.0599613696448501E-006
 uqx(4,3)=   8.6862085239109360E-005
 uqy(4,3)=   2.0599613696448501E-006
 uqx(5,1)=  -1.0375094070201458E-006
 uqy(5,1)=   6.3237707016312353E-007
 uqx(5,2)=   2.1227146117082909E-005
 uqy(5,2)=  -3.9194762329879704E-006
 uqx(5,3)=   3.6782592446672790E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   2.1227146117082909E-005
 uqy(5,4)=   3.9194762329879704E-006
 endif
 if (ntheta.eq.          32 )then
 uqx(2,1)=  -2.7165053726118125E-002
 uqy(2,1)=   1.5795785042915584E-002
 uqx(3,1)=  -1.6284072821584701E-003
 uqy(3,1)=   9.4937741985855911E-004
 uqx(3,2)=  -4.1365569658291003E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -9.4275086002081651E-005
 uqy(4,1)=   5.5133949583323129E-005
 uqx(4,2)=   8.3124866104339803E-005
 uqy(4,2)=   5.5418256663322732E-007
 uqx(4,3)=   8.3124866104339803E-005
 uqy(4,3)=  -5.5418256663322732E-007
 uqx(5,1)=  -2.2521020897502563E-006
 uqy(5,1)=   1.3448724048418223E-006
 uqx(5,2)=   2.3057341827919729E-005
 uqy(5,2)=  -4.0328689250170243E-006
 uqx(5,3)=   4.0936851974362359E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   2.3057341827919729E-005
 uqy(5,4)=   4.0328689250170243E-006
 endif
 if (ntheta.eq.          33 )then
 uqx(2,1)=  -2.8394507311114123E-002
 uqy(2,1)=   1.6506639994148909E-002
 uqx(3,1)=  -1.8211810559783837E-003
 uqy(3,1)=   1.0614214443389939E-003
 uqx(3,2)=  -5.2048939807506744E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -1.1573717794315219E-004
 uqy(4,1)=   6.7639543136682284E-005
 uqx(4,2)=   7.7444096861185304E-005
 uqy(4,2)=   3.6718395435403279E-006
 uqx(4,3)=   7.7444096861185304E-005
 uqy(4,3)=  -3.6718395435403279E-006
 uqx(5,1)=  -3.8635344591073607E-006
 uqy(5,1)=   2.2886040617191943E-006
 uqx(5,2)=   2.4708866548575281E-005
 uqy(5,2)=  -4.0429171197937705E-006
 uqx(5,3)=   4.5083869993072878E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   2.4708866548575281E-005
 uqy(5,4)=   4.0429171197937705E-006
 endif
 if (ntheta.eq.          34 )then
 uqx(2,1)=  -2.9606119757081575E-002
 uqy(2,1)=   1.7207073771532421E-002
 uqx(3,1)=  -2.0246140995809430E-003
 uqy(3,1)=   1.1796309529101220E-003
 uqx(3,2)=  -6.3189666482728031E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -1.4012011579774998E-004
 uqy(4,1)=   8.1840758239034442E-005
 uqx(4,2)=   6.9733069063882344E-005
 uqy(4,2)=   7.3042142251150912E-006
 uqx(4,3)=   6.9733069063882344E-005
 uqy(4,3)=  -7.3042142251150912E-006
 uqx(5,1)=  -5.9371540976553514E-006
 uqy(5,1)=   3.5014592535807844E-006
 uqx(5,2)=   2.6118840472383628E-005
 uqy(5,2)=  -3.9294115964497450E-006
 uqx(5,3)=   4.9153506206729652E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   2.6118840472383628E-005
 uqy(5,4)=   3.9294115964497450E-006
 endif
 if (ntheta.eq.          35 )then
 uqx(2,1)=  -3.0797523100268738E-002
 uqy(2,1)=   1.7895730940821176E-002
 uqx(3,1)=  -2.2382540294391268E-003
 uqy(3,1)=   1.3037421918120587E-003
 uqx(3,2)=  -7.4711417361237199E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -1.6755595008347896E-004
 uqy(4,1)=   9.7813604281450569E-005
 uqx(4,2)=   5.9932702518829840E-005
 uqy(4,2)=   1.1455790873520258E-005
 uqx(4,3)=   5.9932702518829840E-005
 uqy(4,3)=  -1.1455790873520258E-005
 uqx(5,1)=  -8.5401504834948070E-006
 uqy(5,1)=   5.0223706662768778E-006
 uqx(5,2)=   2.7224523372654747E-005
 uqy(5,2)=  -3.6726462538006687E-006
 uqx(5,3)=   5.3074331286291213E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   2.7224523372654747E-005
 uqy(5,4)=   3.6726462538006687E-006
 endif
 if (ntheta.eq.          36 )then
 uqx(2,1)=  -3.1966716465102969E-002
 uqy(2,1)=   1.8571468101473945E-002
 uqx(3,1)=  -2.4616044026482468E-003
 uqy(3,1)=   1.4334661021206542E-003
 uqx(3,2)=  -8.6537658482697000E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -1.9815653042156122E-004
 uqy(4,1)=   1.1562235872043572E-004
 uqx(4,2)=   4.8011092077000848E-005
 uqy(4,2)=   1.6124570727727658E-005
 uqx(4,3)=   4.8011092077000848E-005
 uqy(4,3)=  -1.6124570727727658E-005
 uqx(5,1)=  -1.1740547404553983E-005
 uqy(5,1)=   6.8907327695770384E-006
 uqx(5,2)=   2.7964750724745879E-005
 uqy(5,2)=  -3.2538530853457498E-006
 uqx(5,3)=   5.6775286209476899E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   2.7964750724745879E-005
 uqy(5,4)=   3.2538530853457498E-006
 endif
 if (ntheta.eq.          37 )then
 uqx(2,1)=  -3.3112050475170129E-002
 uqy(2,1)=   1.9233344691700210E-002
 uqx(3,1)=  -2.6941358058514904E-003
 uqy(3,1)=   1.5684947416610797E-003
 uqx(3,2)=  -9.8592993078926877E-004
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -2.3201321158246297E-004
 uqy(4,1)=   1.3531940849307539E-004
 uqx(4,2)=   3.3962308412436504E-005
 uqy(4,2)=   2.1302471188196761E-005
 uqx(4,3)=   3.3962308412436504E-005
 uqy(4,3)=  -2.1302471188196761E-005
 uqx(5,1)=  -1.5606230149174467E-005
 uqy(5,1)=   9.1458389728534019E-006
 uqx(5,2)=   2.8281220212318946E-005
 uqy(5,2)=  -2.6555809419398820E-006
 uqx(5,3)=   6.0187182421142921E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   2.8281220212318946E-005
 uqy(5,4)=   2.6555809419398820E-006
 endif
 if (ntheta.eq.          38 )then
 uqx(2,1)=  -3.4232205798660847E-002
 uqy(2,1)=   1.9880610434395575E-002
 uqx(3,1)=  -2.9352961368716711E-003
 uqy(3,1)=   1.7085072338070947E-003
 uqx(3,2)=  -1.1080428802065397E-003
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -2.6919708406927828E-004
 uqy(4,1)=   1.5694539438150613E-004
 uqx(4,2)=   1.7804615026431277E-005
 uqy(4,2)=   2.6975847572205446E-005
 uqx(4,3)=   1.7804615026431277E-005
 uqy(4,3)=  -2.6975847572205446E-005
 uqx(5,1)=  -2.0204041783697062E-005
 uqy(5,1)=   1.1826359386599187E-005
 uqx(5,2)=   2.8119591278941403E-005
 uqy(5,2)=  -1.8620086446919627E-006
 uqx(5,3)=   6.3243999065080845E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   2.8119591278941403E-005
 uqy(5,4)=   1.8620086446919627E-006
 endif
 if (ntheta.eq.          39 )then
 uqx(2,1)=  -3.5326167618101825E-002
 uqy(2,1)=   2.0512690458548651E-002
 uqx(3,1)=  -3.1845198599577628E-003
 uqy(3,1)=   1.8531751175220020E-003
 uqx(3,2)=  -1.2310158425249287E-003
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -3.0975963923662820E-004
 uqy(4,1)=   1.8052960590128392E-004
 uqx(4,2)=  -4.2173564811337383E-007
 uqy(4,2)=   3.3126099006390713E-005
 uqx(4,3)=  -4.2173564811337383E-007
 uqy(4,3)=  -3.3126099006390713E-005
 uqx(5,1)=  -2.5598973845110378E-005
 uqy(5,1)=   1.4969873731263042E-005
 uqx(5,2)=   2.7430377119402739E-005
 uqy(5,2)=  -8.5918819943473313E-007
 uqx(5,3)=   6.5883953052368218E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   2.7430377119402739E-005
 uqy(5,4)=   8.5918819943473313E-007
 endif
 if (ntheta.eq.          40 )then
 uqx(2,1)=  -3.6393197634710793E-002
 uqy(2,1)=   2.1129169025429276E-002
 uqx(3,1)=  -3.4412361121123137E-003
 uqy(3,1)=   2.0021670288804124E-003
 uqx(3,2)=  -1.3541879904438935E-003
 uqy(3,2)=   0.0000000000000000     
 uqx(4,1)=  -3.5373377446867312E-004
 uqy(4,1)=   2.0609057181782231E-004
 uqx(4,2)=  -2.0656975024795037E-005
 uqy(4,2)=   3.9730323287137485E-005
 uqx(4,3)=  -2.0656975024795037E-005
 uqy(4,3)=  -3.9730323287137485E-005
 uqx(5,1)=  -3.1853467769817504E-005
 uqy(5,1)=   1.8612468714184206E-005
 uqx(5,2)=   2.6169623236922224E-005
 uqy(5,2)=   3.6478171645793533E-007
 uqx(5,3)=   6.8050334755084737E-005
 uqy(5,3)=   0.0000000000000000     
 uqx(5,4)=   2.6169623236922224E-005
 uqy(5,4)=  -3.6478171645793533E-007
 endif


  !Update coordinates
  do i=1,ndim/2

    rx=Coords(i,1)
    ry=Coords(i,2)

    delta(1)=(0.d0,0.d0)
    delta(2)=(0.d0,0.d0)

    do h=1,5,2
      do n=2,5
        do m=1,n-1

          vk(1)=(n-1)*Gn(h,1)+(m-1)*Gn(h+2,1)
          vk(2)=(n-1)*Gn(h,2)+(m-1)*Gn(h+2,2)

          ux=cos((h-1)*pi/3.0_dp)*uqx(n,m)-sin((h-1)*pi/3.0_dp)*uqy(n,m)
          uy=sin((h-1)*pi/3.0_dp)*uqx(n,m)+cos((h-1)*pi/3.0_dp)*uqy(n,m)

          delta(1)=delta(1)-2*ux*sin(vk(1)*rx+vk(2)*ry)
          delta(2)=delta(2)-2*uy*sin(vk(1)*rx+vk(2)*ry)

        end do
      end do
    end do

    Coords(i,1)=rx+delta(1)/2
    Coords(i,2)=ry+delta(2)/2

    Coords(i+ndim/2,1)=-Coords(i,1)!rx-delta(1)/2
    Coords(i+ndim/2,2)=Coords(i,2) !ry-delta(2)/2
   
  enddo

  Nk=0
  do i=ndim/2+1,3*ndim/4
    do j=ndim/2+1,3*ndim/4
      if(sqrt((Coords(i,1) - Coordsnr(j,1))**2 + (Coords(i,2) - Coordsnr(j,2))**2).lt.0.25_dp)then
        Nk = Nk+1
        Coordstemp(j,:) = Coords(i,:)
        exit
      endif 
    enddo
  enddo
  
  Nk=Nk+1
  Coordstemp(3*ndim/4,1) = Coords(3*ndim/4,1)
  Coordstemp(3*ndim/4,2) = -Coords(3*ndim/4,2)
  Coordstemp(3*ndim/4,3) = Coords(3*ndim/4,3)
  
  if(Nk.ne.ndim/4)then
    write(*,*) 'ERROR RELAXATION', Nk
  endif

  Nk=0
  do i=3*ndim/4+1,ndim
    do j=3*ndim/4+1,ndim
      if(sqrt((Coords(i,1) - Coordsnr(j,1))**2 + (Coords(i,2) - Coordsnr(j,2))**2).lt.0.25_dp)then
        Nk = Nk+1
        Coordstemp(j,:) = Coords(i,:)
        exit
      endif 
    enddo
  enddo
  
  Nk=Nk+1
  Coordstemp(ndim,1) = Coords(ndim,1)
  Coordstemp(ndim,2) = -Coords(ndim,2)
  Coordstemp(ndim,3) = Coords(ndim,3)
  
  if(Nk.ne.ndim/4)then
    write(*,*) 'ERROR RELAXATION', Nk
  endif

  Coords(ndim/2+1:ndim,:) = Coordstemp(ndim/2+1:ndim,:)


end subroutine LatticeRelaxationCarr  

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


! Postprocessing compatibility namespace (shares the dev physical settings).
! Compatibility routines grouped by the same responsibilities as FORGE dev.
module ResponseGeometry
use ResponseSetup
use lapack_routines
implicit none
contains

subroutine WignerSeitz(coord,t1,t2,t3,cs,sn)

 real(dp), intent(in)    :: t1(2), t2(2), t3(2), cs, sn
 real(dp), intent(inout) :: coord(ndim,3)

 integer(dp) :: na1, nb1, na2, nb2, n1, n2, nind, nrad
 real(dp)    :: rMax, rTemp1, rTemp2, rTemp3, an1, an2, bn1, bn2

nrad=3*ntheta

! initialize running parameters

nind=0
na1=0
nb1=0
na2=0
nb2=0

 rMax = 3.0_dp*ntheta**2+3.0_dp*ntheta+1.0_dp

 do n1=-nrad,nrad
   do n2=-nrad,nrad

! Find coordinate of site A layer 1 is within 1st Wigner-Seitz cell
! n1*a1+n2*a2

     rTemp1=abs(n1*(3.0_dp*ntheta+1)+n2*(3.0_dp*ntheta+2))+ 1e-6
     rTemp2=abs(n1-n2*(3.0_dp*ntheta+1))+ 1e-6
     rTemp3=abs(n1*(3.0_dp*ntheta+2)+n2)+ 1e-6

     if(rTemp1 < rMax .and. rTemp2 < rMax .and. rTemp3 < rMax)then
       na1=na1+1
       nind=nind+1
       coord(nind,1)=n1*a1(1)+n2*a2(1)
       coord(nind,2)=n1*a1(2)+n2*a2(2)
       coord(nind,3)=0.0_dp
     end if

   end do
 end do

 do n1=-nrad,nrad
   do n2=-nrad,nrad

! Find coordinate of site B layer 1 is within 1st Wigner-Seitz cell
! n1*a1+n2*a2

     bn1=n1+1.0_dp/3.0_dp
     bn2=n2+1.0_dp/3.0_dp

     rTemp1=abs(bn1*(3.0_dp*ntheta+1)+bn2*(3.0_dp*ntheta+2))+ 1e-6
     rTemp2=abs(bn1-bn2*(3.0_dp*ntheta+1))+ 1e-6
     rTemp3=abs(bn1*(3.0_dp*ntheta+2)+bn2)+ 1e-6

     if(rTemp1 < rMax .and. rTemp2 < rMax .and. rTemp3 < rMax)then
       nb1=nb1+1
       nind=nind+1
       coord(nind,1)=bn1*a1(1)+bn2*a2(1)
       coord(nind,2)=bn1*a1(2)+bn2*a2(2)
       coord(nind,3)=0.0_dp
     end if

   end do
 end do

!!! Include B atom at zone boundary

nind=nind+1

coord(nind,1)=(t2(1)+t3(1))/3.0_dp
coord(nind,2)=(t2(2)+t3(2))/3.0_dp
coord(nind,3)=0.0_dp

 do n1=-nrad,nrad
   do n2=-nrad,nrad

! Find coordinate of site A layer 2 is within 1st Wigner-Seitz cell
! n1*a1+n2*a2

     an1=n1*(cs-sn/sq3)-2*n2*sn/sq3
     an2=n2*(cs+sn/sq3)+2*n1*sn/sq3

     rTemp1=abs(an1*(3.0_dp*ntheta+1)+an2*(3.0_dp*ntheta+2))+ 1e-6
     rTemp2=abs(an1-an2*(3.0_dp*ntheta+1))+ 1e-6
     rTemp3=abs(an1*(3.0_dp*ntheta+2)+an2)+ 1e-6

     if(rTemp1 < rMax .and. rTemp2 .LT. rMax .and. rTemp3 < rMax)then
       na2=na2+1
       nind=nind+1
       coord(nind,1)=an1*a1(1)+an2*a2(1)
       coord(nind,2)=an1*a1(2)+an2*a2(2)
       coord(nind,3)=tz
     end if

   end do
 end do

! Find coordinate of site B layer 2 is within 1st Wigner-Seitz cell
! n1*a1+n2*a2

 do n1=-nrad,nrad
   do n2=-nrad,nrad

     bn1=(n1+1.0_dp/3.0_dp)*(cs-sn/sq3)-2*(n2+1.0_dp/3.0_dp)*sn/sq3
     bn2=(n2+1.0_dp/3.0_dp)*(cs+sn/sq3)+2*(n1+1.0_dp/3.0_dp)*sn/sq3

     rTemp1=abs(bn1*(3.0_dp*ntheta+1)+bn2*(3.0_dp*ntheta+2))+ 1e-6
     rTemp2=abs(bn1-bn2*(3.0_dp*ntheta+1))+ 1e-6
     rTemp3=abs(bn1*(3.0_dp*ntheta+2)+bn2)+ 1e-6

     if(rTemp1 < rMax .AND. rTemp2 < rMax .AND. rTemp3 < rMax)then
       nb2=nb2+1
       nind=nind+1
       coord(nind,1)=bn1*a1(1)+bn2*a2(1)
       coord(nind,2)=bn1*a1(2)+bn2*a2(2)
       coord(nind,3)=tz
     end if

   end do
 end do

!!! Include B atom at zone boundary

 nind=nind+1

 coord(nind,1)=(t1(1)+t2(1))/3.0_dp
 coord(nind,2)=(t1(2)+t2(2))/3.0_dp
 coord(nind,3)=tz

end subroutine WignerSeitz





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine WignerSeitz6(coord,t1,t2,t3,cs,sn)


  real(dp), intent(in)    :: t1(2), t2(2), t3(2), cs, sn

  real(dp), intent(inout) :: coord(ndim,3)



  integer(dp) :: na1, nb1, na2, nb2, n1, n2, nind, nrad

  real(dp)    :: rMax, rTemp1, rTemp2, rTemp3, an1, an2, bn1, bn2



 nrad=3*ntheta



 ! initialize running parameters



 nind=0

 na1=0

 nb1=0

 na2=0

 nb2=0



  rMax = 3.0_dp*ntheta**2+3.0_dp*ntheta+1.0_dp



  do n1=-nrad,nrad

    do n2=-nrad,nrad



 ! Find coordinate of site A layer 1 is within 1st Wigner-Seitz cell

 ! n1*a1+n2*a2

      an1=(n1+1.0_dp/3.0_dp)

      an2=(n2-2.0_dp/3.0_dp)



      rTemp1=abs(an1*(3.0_dp*ntheta+1)+an2*(3.0_dp*ntheta+2))+ 1e-6

      rTemp2=abs(an1-an2*(3.0_dp*ntheta+1))+ 1e-6

      rTemp3=abs(an1*(3.0_dp*ntheta+2)+an2)+ 1e-6



      if(rTemp1 < rMax .and. rTemp2 < rMax .and. rTemp3 < rMax)then

        na1=na1+1

        nind=nind+1

        coord(nind,1)=(n1+1.0_dp/3.0_dp)*a1(1)+(n2-2.0_dp/3.0_dp)*a2(1)

        coord(nind,2)=(n1+1.0_dp/3.0_dp)*a1(2)+(n2-2.0_dp/3.0_dp)*a2(2)

        coord(nind,3)=0.0_dp

      end if



    end do

  end do



  !! Include A atom at zone boundary



 nind=nind+1



 coord(nind,1)=(-t2(1)-t1(1))/3.0_dp

 coord(nind,2)=(-t2(2)-t1(2))/3.0_dp

 coord(nind,3)=0.0_dp



  do n1=-nrad,nrad

    do n2=-nrad,nrad



 ! Find coordinate of site B layer 1 is within 1st Wigner-Seitz cell

 ! n1*a1+n2*a2



      bn1=n1+2.0_dp/3.0_dp

      bn2=n2-1.0_dp/3.0_dp



      rTemp1=abs(bn1*(3.0_dp*ntheta+1)+bn2*(3.0_dp*ntheta+2))+ 1e-6

      rTemp2=abs(bn1-bn2*(3.0_dp*ntheta+1))+ 1e-6

      rTemp3=abs(bn1*(3.0_dp*ntheta+2)+bn2)+ 1e-6



      if(rTemp1 < rMax .and. rTemp2 < rMax .and. rTemp3 < rMax)then

        nb1=nb1+1

        nind=nind+1

        coord(nind,1)=bn1*a1(1)+bn2*a2(1)

        coord(nind,2)=bn1*a1(2)+bn2*a2(2)

        coord(nind,3)=0.0_dp

      end if



    end do

  end do



 !! Include B atom at zone boundary



 nind=nind+1



 coord(nind,1)=(t2(1)+t1(1))/3.0_dp

 coord(nind,2)=(t2(2)+t1(2))/3.0_dp

 coord(nind,3)=0.0_dp



  do n1=-nrad,nrad

    do n2=-nrad,nrad



 ! Find coordinate of site A layer 2 is within 1st Wigner-Seitz cell

 ! n1*a1+n2*a2



      an1=(n1+1.0_dp/3.0_dp)*(cs-sn/sq3)-2*(n2-2.0_dp/3.0_dp)*sn/sq3

      an2=(n2-2.0_dp/3.0_dp)*(cs+sn/sq3)+2*(n1+1.0_dp/3.0_dp)*sn/sq3



      rTemp1=abs(an1*(3.0_dp*ntheta+1)+an2*(3.0_dp*ntheta+2))+ 1e-6

      rTemp2=abs(an1-an2*(3.0_dp*ntheta+1))+ 1e-6

      rTemp3=abs(an1*(3.0_dp*ntheta+2)+an2)+ 1e-6



      if(rTemp1 < rMax .and. rTemp2 .LT. rMax .and. rTemp3 < rMax)then

        na2=na2+1

        nind=nind+1

        coord(nind,1)=an1*a1(1)+an2*a2(1)

        coord(nind,2)=an1*a1(2)+an2*a2(2)

        coord(nind,3)=tz

      end if



    end do

  end do



  !! Include A atom at zone boundary



 nind=nind+1



 coord(nind,1)=(t1(1)+t2(1))/3.0_dp

 coord(nind,2)=(t1(2)+t2(2))/3.0_dp

 coord(nind,3)=tz



 ! Find coordinate of site B layer 2 is within 1st Wigner-Seitz cell

 ! n1*a1+n2*a2



  do n1=-nrad,nrad

    do n2=-nrad,nrad



      bn1=(n1+2.0_dp/3.0_dp)*(cs-sn/sq3)-2*(n2-1.0_dp/3.0_dp)*sn/sq3

      bn2=(n2-1.0_dp/3.0_dp)*(cs+sn/sq3)+2*(n1+2.0_dp/3.0_dp)*sn/sq3



      rTemp1=abs(bn1*(3.0_dp*ntheta+1)+bn2*(3.0_dp*ntheta+2))+ 1e-6

      rTemp2=abs(bn1-bn2*(3.0_dp*ntheta+1))+ 1e-6

      rTemp3=abs(bn1*(3.0_dp*ntheta+2)+bn2)+ 1e-6



      if(rTemp1 < rMax .AND. rTemp2 < rMax .AND. rTemp3 < rMax)then

        nb2=nb2+1

        nind=nind+1

        coord(nind,1)=bn1*a1(1)+bn2*a2(1)

        coord(nind,2)=bn1*a1(2)+bn2*a2(2)

        coord(nind,3)=tz

      end if



    end do

  end do



 !!! Include B atom at zone boundary



  nind=nind+1



  coord(nind,1)=(-t1(1)-t2(1))/3.0_dp

  coord(nind,2)=(-t1(2)-t2(2))/3.0_dp

  coord(nind,3)=tz



end subroutine WignerSeitz6



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Relaxation_old(coord,Gn)

 real(dp), intent(inout) :: coord(ndim,3)
 real(dp), intent(in)    :: Gn(8,2)

 integer(dp) :: i, m, n, ncount, nn
 real(dp) :: ux, uy, rx, ry
 real(dp) :: uqx(5,5), uqy(5,5), delta(2)
 real(dp) :: vk(2)
 real(dp), allocatable :: uq(:,:)

 allocate(uq(ndim,2))

uqx=0.0_dp
uqy=0.0_dp

if(ntheta.EQ.2)then
uqx(2,1)=-0.0003688729746213733
uqy(2,1)=0.00023674102395776823
uqx(3,1)=-2.930162659354381e-07
uqy(3,1)=1.880805670248788e-07
uqx(3,2)=-4.0348301280806025e-07
uqy(3,2)=0.0_dp
uqx(4,1)=-3.103595720911555e-10
uqy(4,1)=1.9921496057160015e-10
uqx(4,2)=-5.555214771937626e-10
uqy(4,2)=1.2701266883624605e-10
uqx(4,3)=-5.555214771937626e-10
uqy(4,3)=-1.2701266883624605e-10
uqx(5,1)=-3.6798131426275956e-13
uqy(5,1)=2.3619762123031145e-13
uqx(5,2)=-8.303523356604859e-13
uqy(5,2)=2.7790899851632874e-13
uqx(5,3)=-9.564169350311786e-13
uqy(5,3)=0.0_dp
uqx(5,4)=-8.303523356604859e-13
uqy(5,4)=-2.7790899851632874e-13
endif
if(ntheta.EQ.3)then
uqx(2,1)=-0.0017606949668591743
uqy(2,1)=0.001096087467910831
uqx(3,1)=-6.665490380552498e-06
uqy(3,1)=4.15132641511213e-06
uqx(3,2)=-8.962118426312275e-06
uqy(3,2)=0.0_dp
uqx(4,1)=-3.374281857549946e-08
uqy(4,1)=2.1016000503612443e-08
uqx(4,2)=-5.8677599481932206e-08
uqy(4,2)=1.3531357298462511e-08
uqx(4,3)=-5.8677599481932206e-08
uqy(4,3)=-1.3531357298462511e-08
uqx(5,1)=-1.9277559881782495e-10
uqy(5,1)=1.2005800437758387e-10
uqx(5,2)=-4.2095892504204485e-10
uqy(5,2)=1.4171817228651717e-10
uqx(5,3)=-4.79235073634704e-10
uqy(5,3)=0.0_dp
uqx(5,4)=-4.2095892504204485e-10
uqy(5,4)=-1.4171817228651717e-10
endif
if(ntheta.EQ.4)then
uqx(2,1)=-0.002902656470153326
uqy(2,1)=0.001776536901182341
uqx(3,1)=-1.8131963434478182e-05
uqy(3,1)=1.1103919737259968e-05
uqx(3,2)=-2.399383909925256e-05
uqy(3,2)=0.0_dp
uqx(4,1)=-1.5171149243510273e-07
uqy(4,1)=9.291155272630708e-08
uqx(4,2)=-2.586971229650243e-07
uqy(4,2)=5.993758537213166e-08
uqx(4,3)=-2.586971229650243e-07
uqy(4,3)=-5.993758537213166e-08
uqx(5,1)=-1.4355267714309013e-09
uqy(5,1)=8.790753953291303e-10
uqx(5,2)=-3.0631664980077603e-09
uqy(5,2)=1.0346849951479528e-09
uqx(5,3)=-3.462319317853007e-09
uqy(5,3)=0.0_dp
uqx(5,4)=-3.0631664980077603e-09
uqy(5,4)=-1.0346849951479528e-09
endif
if(ntheta.EQ.5)then
uqx(2,1)=-0.004311382273032596
uqy(2,1)=0.0026099257526777447
uqx(3,1)=-4.015514370958257e-05
uqy(3,1)=2.4325749055854244e-05
uqx(3,2)=-5.229034708204007e-05
uqy(3,2)=0.0_dp
uqx(4,1)=-5.018320220897382e-07
uqy(4,1)=3.040253371439883e-07
uqx(4,2)=-8.386608841906334e-07
uqy(4,2)=1.950837119571707e-07
uqx(4,3)=-8.386608841906334e-07
uqy(4,3)=-1.950837119571707e-07
uqx(5,1)=-7.108836115162087e-09
uqy(5,1)=4.306334963933745e-09
uqx(5,2)=-1.4805138782519942e-08
uqy(5,2)=5.013360981832663e-09
uqx(5,3)=-1.6611640114265272e-08
uqy(5,3)=0.0_dp
uqx(5,4)=-1.4805138782519942e-08
uqy(5,4)=-5.013360981832663e-09
endif
if(ntheta.EQ.6)then
uqx(2,1)=-0.005969705915838723
uqy(2,1)=0.003586090624014734
uqx(3,1)=-7.745829739097424e-05
uqy(3,1)=4.657056067287355e-05
uqx(3,2)=-9.913176821488675e-05
uqy(3,2)=0.0_dp
uqx(4,1)=-1.350887244346724e-06
uqy(4,1)=8.122623735641327e-07
uqx(4,2)=-2.2090102714001943e-06
uqy(4,2)=5.158716841134904e-07
uqx(4,3)=-2.2090102714001943e-06
uqy(4,3)=-5.158716841134904e-07
uqx(5,1)=-2.6771652167817564e-08
uqy(5,1)=1.6095542724748687e-08
uqx(5,2)=-5.430792193700543e-08
uqy(5,2)=1.8431133037986742e-08
uqx(5,3)=-6.044789289403583e-08
uqy(5,3)=0.0_dp
uqx(5,4)=-5.430792193700543e-08
uqy(5,4)=-1.8431133037986742e-08
endif
if(ntheta.EQ.7)then
uqx(2,1)=-0.007854823723483664
uqy(2,1)=0.0046916376939870035
uqx(3,1)=-0.00013519946472068375
uqy(3,1)=8.083528462727111e-05
uqx(3,2)=-0.00016975181256366468
uqy(3,2)=0.0_dp
uqx(4,1)=-3.1324906494385004e-06
uqy(4,1)=1.8730846215242352e-06
uqx(4,2)=-5.002321889185792e-06
uqy(4,2)=1.1730320003764776e-06
uqx(4,3)=-5.002321889185792e-06
uqy(4,3)=-1.1730320003764776e-06
uqx(5,1)=-8.268500554445054e-08
uqy(5,1)=4.943605235578102e-08
uqx(5,2)=-1.6303025257307402e-07
uqy(5,2)=5.545235390537088e-08
uqx(5,3)=-1.7987656889364531e-07
uqy(5,3)=0.0_dp
uqx(5,4)=-1.6303025257307402e-07
uqy(5,4)=-5.545235390537088e-08
endif
if(ntheta.EQ.8)then
uqx(2,1)=-0.00993810491283713
uqy(2,1)=0.0059098453909551965
uqx(3,1)=-0.0002186142755737399
uqy(3,1)=0.00013015245894254537
uqx(3,2)=-0.00026875584006595866
uqy(3,2)=0.0_dp
uqx(4,1)=-6.480157897394943e-06
uqy(4,1)=3.85841764396702e-06
uqx(4,2)=-1.0085970871722486e-05
uqy(4,2)=2.3755731921643997e-06
uqx(4,3)=-1.0085970871722486e-05
uqy(4,3)=-2.3755731921643997e-06
uqx(5,1)=-2.1939330558573372e-07
uqy(5,1)=1.3061506709174655e-07
uqx(5,2)=-4.1965184358195726e-07
uqy(5,2)=1.430645718761784e-07
uqx(5,3)=-4.586390390493329e-07
uqy(5,3)=0.0_dp
uqx(5,4)=-4.1965184358195726e-07
uqy(5,4)=-1.430645718761784e-07
endif
if(ntheta.EQ.9)then
uqx(2,1)=-0.012185588828143494
uqy(2,1)=0.007220966695102083
uqx(3,1)=-0.00033260151238692366
uqy(3,1)=0.00019735023443253702
uqx(3,2)=-0.0003995486377704783
uqy(3,2)=0.0_dp
uqx(4,1)=-1.2238548704210078e-05
uqy(4,1)=7.262751974716624e-06
uqx(4,2)=-1.8531988726911148e-05
uqy(4,2)=4.385357805281378e-06
uqx(4,3)=-1.8531988726911148e-05
uqy(4,3)=-4.385357805281378e-06
uqx(5,1)=-5.156154951107497e-07
uqy(5,1)=3.0594317571121564e-07
uqx(5,2)=-9.55294885896738e-07
uqy(5,2)=3.264465793803411e-07
uqx(5,3)=-1.0335413851090244e-06
uqy(5,3)=0.0_dp
uqx(5,4)=-9.55294885896738e-07
uqy(5,4)=-3.264465793803411e-07
endif
if(ntheta.EQ.10)then
uqx(2,1)=-0.014559262736192947
uqy(2,1)=0.008602980245994919
uqx(3,1)=-0.00048130620605092467
uqy(3,1)=0.0002848118672289576
uqx(3,2)=-0.0005638718113666972
uqy(3,2)=0.0_dp
uqx(4,1)=-2.1443881780707192e-05
uqy(4,1)=1.269128122010679e-05
uqx(4,2)=-3.1540184842605977e-05
uqy(4,2)=7.50038681991021e-06
uqx(4,3)=-3.1540184842605977e-05
uqy(4,3)=-7.50038681991021e-06
uqx(5,1)=-1.0963582385959045e-06
uqy(5,1)=6.487770601854134e-07
uqx(5,2)=-1.9651901070438073e-06
uqy(5,2)=6.732111753546448e-07
uqx(5,3)=-2.103770247456669e-06
uqy(5,3)=0.0_dp
uqx(5,4)=-1.9651901070438073e-06
uqy(5,4)=-6.732111753546448e-07
endif
if(ntheta.EQ.11)then
uqx(2,1)=-0.017019014749540734
uqy(2,1)=0.010032725521765675
uqx(3,1)=-0.0006677701296776364
uqy(3,1)=0.0003942751299028151
uqx(3,2)=-0.000761545097275991
uqy(3,2)=0.0_dp
uqx(4,1)=-3.5271199903725905e-05
uqy(4,1)=2.082882391798765e-05
uqx(4,2)=-5.032398576258808e-05
uqy(4,2)=1.202842921648509e-05
uqx(4,3)=-5.032398576258808e-05
uqy(4,3)=-1.202842921648509e-05
uqx(5,1)=-2.1421421129757634e-06
uqy(5,1)=1.2648269448737656e-06
uqx(5,2)=-3.7120170614227704e-06
uqy(5,2)=1.2748598470229163e-06
uqx(5,3)=-3.930711929429223e-06
uqy(5,3)=0.0_dp
uqx(5,4)=-3.7120170614227704e-06
uqy(5,4)=-1.2748598470229163e-06
endif
if(ntheta.EQ.12)then
uqx(2,1)=-0.019524958181451703
uqy(2,1)=0.011487244375642828
uqx(3,1)=-0.0008937126804577912
uqy(3,1)=0.0005267078513122274
uqx(3,2)=-0.0009904654389958113
uqy(3,2)=0.0_dp
uqx(4,1)=-5.495527348796896e-05
uqy(4,1)=3.239356777699126e-05
uqx(4,2)=-7.598082680854464e-05
uqy(4,2)=1.825569088821883e-05
uqx(4,3)=-7.598082680854464e-05
uqy(4,3)=-1.825569088821883e-05
uqx(5,1)=-3.892204427996963e-06
uqy(5,1)=2.293932119077245e-06
uqx(5,2)=-6.517721327411936e-06
uqy(5,2)=2.244282599798575e-06
uqx(5,3)=-6.825897690417241e-06
uqy(5,3)=0.0_dp
uqx(5,4)=-6.517721327411936e-06
uqy(5,4)=-2.244282599798575e-06
endif
if(ntheta.EQ.13)then
uqx(2,1)=-0.02203971289161037
uqy(2,1)=0.012945089479999937
uqx(3,1)=-0.001159475806587175
uqy(3,1)=0.0006822786428372063
uqx(3,2)=-0.0012468548139774377
uqy(3,2)=0.0_dp
uqx(4,1)=-8.169970517655794e-05
uqy(4,1)=4.8084156180646557e-05
uqx(4,2)=-0.0001093744119432367
uqy(4,2)=2.6416874142175212e-05
uqx(4,3)=-0.0001093744119432367
uqy(4,3)=-2.6416874142175212e-05
uqx(5,1)=-6.639802204729499e-06
uqy(5,1)=3.907231157592452e-06
uqx(5,2)=-1.0743817669032263e-05
uqy(5,2)=3.7092121816072646e-06
uqx(5,3)=-1.1128290721222887e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-1.0743817669032263e-05
uqy(5,4)=-3.7092121816072646e-06
endif
if(ntheta.EQ.14)then
uqx(2,1)=-0.024530256594078644
uqy(2,1)=0.014387378212308049
uqx(3,1)=-0.0014641263855930343
uqy(3,1)=0.0008604185419698234
uqx(3,2)=-0.0015256885873251322
uqy(3,2)=0.0_dp
uqx(4,1)=-0.00011659167117290876
uqy(4,1)=6.853025116032777e-05
uqx(4,2)=-0.00015105198106743665
uqy(4,2)=3.667243266593343e-05
uqx(4,3)=-0.00015105198106743665
uqy(4,3)=-3.667243266593343e-05
uqx(5,1)=-1.0719859710201104e-05
uqy(5,1)=6.299891069023662e-06
uqx(5,2)=-1.676363535134788e-05
uqy(5,2)=5.802737514994693e-06
uqx(5,3)=-1.7175456228934552e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-1.676363535134788e-05
uqy(5,4)=-5.802737514994693e-06
endif
if(ntheta.EQ.15)then
uqx(2,1)=-0.026969104893425525
uqy(2,1)=0.015798455965871965
uqx(3,1)=-0.001805675239338205
uqy(3,1)=0.0010599496130215184
uqx(3,2)=-0.0018212044309949956
uqy(3,2)=0.0_dp
uqx(4,1)=-0.00016053668592715893
uqy(4,1)=9.42548851073656e-05
uqx(4,2)=-0.00020120823453069478
uqy(4,2)=4.909642780579341e-05
uqx(4,3)=-0.00020120823453069478
uqy(4,3)=-4.909642780579341e-05
uqx(5,1)=-1.6491243899095915e-05
uqy(5,1)=9.68074128635529e-06
uqx(5,2)=-2.4932073125853637e-05
uqy(5,2)=8.652751748689217e-06
uqx(5,3)=-2.5274679964717916e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-2.4932073125853637e-05
uqy(5,4)=-8.652751748689217e-06
endif
if(ntheta.EQ.16)then
uqx(2,1)=-0.02933476915973707
uqy(2,1)=0.017166142817697613
uqx(3,1)=-0.0021813557685715526
uqy(3,1)=0.0012792473587272606
uqx(3,2)=-0.0021273952549505904
uqy(3,2)=0.0_dp
uqx(4,1)=-0.00021422081310178742
uqy(4,1)=0.0001256528407670863
uqx(4,2)=-0.00025969444982031607
uqy(4,2)=6.367535937378761e-05
uqx(4,3)=-0.00025969444982031607
uqy(4,3)=-6.367535937378761e-05
uqx(5,1)=-2.431705664771852e-05
uqy(5,1)=1.4260791613716685e-05
uqx(5,2)=-3.555845401412153e-05
uqy(5,2)=1.2372275334403248e-05
uqx(5,3)=-3.567937813962246e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-3.555845401412153e-05
uqy(5,4)=-1.2372275334403248e-05
endif
if(ntheta.EQ.17)then
uqx(2,1)=-0.03161159838733411
uqy(2,1)=0.018481627262408935
uqx(3,1)=-0.002587908376466652
uqy(3,1)=0.0015164058211511447
uqx(3,2)=-0.0024384184041005663
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0002781002756859611
uqy(4,1)=0.0001629849598205697
uqx(4,2)=-0.0003260618186293725
uqy(4,2)=8.031593146865948e-05
uqx(4,3)=-0.0003260618186293725
uqy(4,3)=-8.031593146865948e-05
uqx(5,1)=-3.454623622703092e-05
uqy(5,1)=2.024254671882617e-05
uqx(5,2)=-4.888639514695799e-05
uqy(5,2)=1.7052063096837917e-05
uqx(5,3)=-4.8573910203114055e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-4.888639514695799e-05
uqy(5,4)=-1.7052063096837917e-05
endif
if(ntheta.EQ.18)then
uqx(2,1)=-0.033789193644226535
uqy(2,1)=0.019739117090906175
uqx(3,1)=-0.003021832398265609
uqy(3,1)=0.0017693833619352496
uqx(3,2)=-0.0027488901953418
uqy(3,2)=0.0_dp
uqx(4,1)=-0.00035241311189526084
uqy(4,1)=0.0002063852171666288
uqx(4,2)=-0.0003996243011955691
uqy(4,2)=9.885854958894279e-05
uqx(4,3)=-0.0003996243011955691
uqy(4,3)=-9.885854958894279e-05
uqx(5,1)=-4.7498791109832806e-05
uqy(5,1)=2.7811458869794964e-05
uqx(5,2)=-6.508227991823797e-05
uqy(5,2)=2.2756121742700963e-05
uqx(5,3)=-6.406738379493807e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-6.508227991823797e-05
uqy(5,4)=-2.2756121742700963e-05
endif
if(ntheta.EQ.19)then
uqx(2,1)=-0.035861593572565124
uqy(2,1)=0.020935362099899803
uqx(3,1)=-0.003479585983636394
uqy(3,1)=0.0020361179886394156
uqx(3,2)=-0.003054065572278459
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0004372052794297127
uqy(4,1)=0.0002558761197913376
uqx(4,2)=-0.00047952745493926924
uqy(4,2)=0.00011909334780127935
uqx(4,3)=-0.00047952745493926924
uqy(4,3)=-0.00011909334780127935
uqx(5,1)=-6.34557265075311e-05
uqy(5,1)=3.7130123742686886e-05
uqx(5,2)=-8.423190528103735e-05
uqy(5,2)=2.9520064581031473e-05
uqx(5,3)=-8.21951573754533e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-8.423190528103735e-05
uqy(5,4)=-2.9520064581031473e-05
endif
if(ntheta.EQ.20)then
uqx(2,1)=-0.03782639325500282
uqy(2,1)=0.02206914176944493
uqx(3,1)=-0.003957729814560495
uqy(3,1)=0.0023146100099506153
uqx(3,2)=-0.0033499211638098473
uqy(3,2)=0.0_dp
uqx(4,1)=-0.000532363868208654
uqy(4,1)=0.00031138816043362787
uqx(4,2)=-0.0005648134487120218
uqy(4,2)=0.00014077624693736233
uqx(4,3)=-0.0005648134487120218
uqy(4,3)=-0.00014077624693736233
uqx(5,1)=-8.26536411922198e-05
uqy(5,1)=4.833519554970959e-05
uqx(5,2)=-0.00010634365464369339
uqy(5,2)=3.7351788455729954e-05
uqx(5,3)=-0.00010292582645481166
uqy(5,3)=0.0_dp
uqx(5,4)=-0.00010634365464369339
uqy(5,4)=-3.7351788455729954e-05
endif
if(ntheta.EQ.21)then
uqx(2,1)=-0.039683905988072816
uqy(2,1)=0.02314078029726673
uqx(3,1)=-0.004453020231804285
uqy(3,1)=0.002602975367839706
uqx(3,2)=-0.0036331670187556417
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0006376517061537844
uqy(4,1)=0.0003727800138387269
uqx(4,2)=-0.0006544766965179579
uqy(4,2)=0.0001636434595389909
uqx(4,3)=-0.0006544766965179579
uqy(4,3)=-0.0001636434595389909
uqx(5,1)=-0.00010528328537327498
uqy(5,1)=6.153660137768897e-05
uqx(5,2)=-0.00013135612165942703
uqy(5,2)=4.6233787024697836e-05
uqx(5,3)=-0.00012617137716973832
uqy(5,3)=0.0_dp
uqx(5,4)=-0.00013135612165942703
uqy(5,4)=-4.6233787024697836e-05
endif
if(ntheta.EQ.22)then
uqx(2,1)=-0.04143642815144943
uqy(2,1)=0.024151722886182475
uqx(3,1)=-0.0049624617327091725
uqy(3,1)=0.002899475489422576
uqx(3,2)=-0.003901211536393849
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0007527396276284611
uqy(4,1)=0.0004398573228967515
uqx(4,2)=-0.0007475079925818254
uqy(4,2)=0.00018742366805470892
uqx(4,3)=-0.0007475079925818254
uqy(4,3)=-0.00018742366805470892
uqx(5,1)=-0.0001314910775426136
uqy(5,1)=7.681846766995464e-05
uqx(5,2)=-0.00015914823875347363
uqy(5,2)=5.6126438082737944e-05
uqx(5,3)=-0.00015179857171594715
uqy(5,3)=0.0_dp
uqx(5,4)=-0.00015914823875347363
uqy(5,4)=-5.6126438082737944e-05
endif
if(ntheta.EQ.23)then
uqx(2,1)=-0.04308763035883389
uqy(2,1)=0.02510418596199579
uqx(3,1)=-0.005483329547043717
uqy(3,1)=0.003202529894557779
uqx(3,2)=-0.004152099585825217
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0008772344173302233
uqy(4,1)=0.0005123889356509687
uqx(4,2)=-0.0008429273067893358
uqy(4,2)=0.0002118476820384999
uqx(4,3)=-0.0008429273067893358
uqy(4,3)=-0.0002118476820384999
uqx(5,1)=-0.0001613825687053204
uqy(5,1)=9.424116901736953e-05
uqx(5,2)=-0.00018955036565760028
uqy(5,2)=6.697172673370301e-05
uqx(5,3)=-0.0001796401856940401
uqy(5,3)=0.0_dp
uqx(5,4)=-0.00018955036565760028
uqy(5,4)=-6.697172673370301e-05
endif
if(ntheta.EQ.24)then
uqx(2,1)=-0.04464207458235026
uqy(2,1)=0.026000880622690153
uqx(3,1)=-0.006013171781584232
uqy(3,1)=0.003510717065719575
uqx(3,2)=-0.004384438316998377
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0010107017141430948
uqy(4,1)=0.0005901201910458327
uqx(4,2)=-0.0009398066180065344
uqy(4,2)=0.00023665572578278698
uqx(4,3)=-0.0009398066180065344
uqy(4,3)=-0.00023665572578278698
uqx(5,1)=-0.00019502698873427135
uqy(5,1)=0.00011384399637522258
uqx(5,2)=-0.0002223552583997825
uqy(5,2)=7.869701860174408e-05
uqx(5,3)=-0.00020950525091773434
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0002223552583997825
uqy(5,4)=-7.869701860174408e-05
endif
if(ntheta.EQ.25)then
uqx(2,1)=-0.04610484423986092
uqy(2,1)=0.026844801459071962
uqx(3,1)=-0.0065497986804440664
uqy(3,1)=0.003822767946492368
uqx(3,2)=-0.004597320221149113
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0011526839626769694
uqy(4,1)=0.000672783314919442
uqx(4,2)=-0.0010372846258475965
uqy(4,2)=0.0002616026744987021
uqx(4,3)=-0.0010372846258475965
uqy(4,3)=-0.0002616026744987021
uqx(5,1)=-0.00023246221336892673
uqy(5,1)=0.0001356480613594293
uqx(5,2)=-0.0002573282525611999
uqy(5,2)=9.121863665062806e-05
uqx(5,3)=-0.0002411878797633425
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0002573282525611999
uqy(5,4)=-9.121863665062806e-05
endif
if(ntheta.EQ.26)then
uqx(2,1)=-0.047481258615900926
uqy(2,1)=0.02763906395972363
uqx(3,1)=-0.007091265279716939
uqy(3,1)=0.004137555698215187
uqx(3,2)=-0.004790249012643126
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0013027140647043248
uqy(4,1)=0.0007601053134247232
uqx(4,2)=-0.0011345752539357542
uqy(4,2)=0.0002864616726489501
uqx(4,3)=-0.0011345752539357542
uqy(4,3)=-0.0002864616726489501
uqx(5,1)=-0.0002736997183935738
uqy(5,1)=0.00015965918640975345
uqx(5,2)=-0.0002942163372629156
uqy(5,2)=0.00010444512093628974
uqx(5,3)=-0.000274474549943679
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0002942163372629156
uqy(5,4)=-0.00010444512093628974
endif
if(ntheta.EQ.27)then
uqx(2,1)=-0.048776353070479665
uqy(2,1)=0.028386605209054172
uqx(3,1)=-0.007635870700194479
uqy(3,1)=0.0044540952921595205
uqx(3,2)=-0.00496306544569023
uqy(3,2)=0.0_dp
uqx(4,1)=-0.001460329223855284
uqy(4,1)=0.0008518160008897586
uqx(4,2)=-0.0012309732077084414
uqy(4,2)=0.0003110283819314306
uqx(4,3)=-0.0012309732077084414
uqy(4,3)=-0.0003110283819314306
uqx(5,1)=-0.0003187300552444502
uqy(5,1)=0.00018587109683292205
uqx(5,2)=-0.0003327567283867869
uqy(5,2)=0.00011828060696595211
uqx(5,3)=-0.00030915014796716294
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0003327567283867869
uqy(5,4)=-0.00011828060696595211
endif
if(ntheta.EQ.28)then
uqx(2,1)=-0.04998906281762806
uqy(2,1)=0.029086805073427164
uqx(3,1)=-0.008182485709112529
uqy(3,1)=0.004771735434979995
uqx(3,2)=-0.005115749172624432
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0016251493713437886
uqy(4,1)=0.0009476939242635575
uqx(4,2)=-0.0013258819305531008
uqy(4,2)=0.00033515809254836736
uqx(4,3)=-0.0013258819305531008
uqy(4,3)=-0.00033515809254836736
uqx(5,1)=-0.0003675432812169361
uqy(5,1)=0.00021427736986447406
uqx(5,2)=-0.0003726973347517757
uqy(5,2)=0.00013263695984855765
uqx(5,3)=-0.00034500664330239087
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0003726973347517757
uqy(5,4)=-0.00013263695984855765
endif
if(ntheta.EQ.29)then
uqx(2,1)=-0.051042686286214146
uqy(2,1)=0.02969507739823804
uqx(3,1)=-0.00873460347786944
uqy(3,1)=0.005092530866296964
uqx(3,2)=-0.005246372723420572
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0017977534684612002
uqy(4,1)=0.0010480801465725953
uqx(4,2)=-0.0014191017341478457
uqy(4,2)=0.0003592056219596041
uqx(4,3)=-0.0014191017341478457
uqy(4,3)=-0.0003592056219596041
uqx(5,1)=-0.00042034060028563143
uqy(5,1)=0.00024499533553742523
uqx(5,2)=-0.0004139705188268175
uqy(5,2)=0.00014755781479141694
uqx(5,3)=-0.0003818914850522088
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0004139705188268175
uqy(5,4)=-0.00014755781479141694
endif
if(ntheta.EQ.300)then !Koshino
uqx(2,1)=-3.5955179495505439E-002
uqy(2,1)=2.0906415366297500E-002
uqx(3,1)=-3.4974142803225836E-003
uqy(3,1)=2.0366806792244314E-003
uqx(3,2)=-3.0623567467657864E-003
uqy(3,2)=0.0000000000000000
uqx(4,1)=-4.4036485598978421E-004
uqy(4,1)=2.5646772695932669E-004
uqx(4,2)=-4.8184834784040605E-004
uqy(4,2)=1.1970911652618624E-004
uqx(4,3)=-4.8184834784040605E-004
uqy(4,3)=-1.1970911652618624E-004
uqx(5,1)=-6.4043997656082756E-005
uqy(5,1)=3.7294216529563565E-005
uqx(5,2)=-8.4819889884739453E-005
uqy(5,2)=2.9733239460941308E-005
uqx(5,3)=-8.2733545432358833E-005
uqy(5,3)=0.0000000000000000
uqx(5,4)=-8.4819889884739453E-005
uqy(5,4)=-2.9733239460941308E-005
endif
if(ntheta.EQ.30)then
uqx(2,1)=-0.05116174430212328
uqy(2,1)=0.029760870864795994
uqx(3,1)=-0.009335082847245551
uqy(3,1)=0.005441585368438033
uqx(3,2)=-0.005330415896682406
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0019874789431051708
uqy(4,1)=0.0011584506658465247
uqx(4,2)=-0.0015129554478849553
uqy(4,2)=0.0003879757319977168
uqx(4,3)=-0.0015129554478849553
uqy(4,3)=-0.0003879757319977168
uqx(5,1)=-0.00047948174355565083
uqy(5,1)=0.0002794059394314642
uqx(5,2)=-0.00045817363144676297
uqy(5,2)=0.00016434148788292767
uqx(5,3)=-0.0004200244526097355
uqy(5,3)=0.0_dp
uqx(5,4)=-0.00045817363144676297
uqy(5,4)=-0.00016434148788292767
endif
if(ntheta.EQ.310)then !original
uqx(1+1,0+1)=-0.0372_dp
uqy(1+1,0+1)=0.02148_dp
uqx(2+1,0+1)=-0.003773_dp
uqy(2+1,0+1)=0.002178_dp
uqx(2+1,1+1)=-0.00324_dp
uqy(2+1,1+1)=0._dp
uqx(3+1,0+1)=-0.0004952_dp
uqy(3+1,0+1)=0.0002859_dp
uqx(3+1,1+1)=-0.0005003_dp
uqy(3+1,1+1)=0.0002456_dp
uqx(3+1,2+1)=-0.0005003_dp
uqy(3+1,2+1)=-0.0002456_dp
endif
if(ntheta.EQ.310)then !Koshino
uqx(2,1)=-3.5955179495505439E-002
uqy(2,1)=2.0906415366297500E-002
uqx(3,1)=-3.4974142803225836E-003
uqy(3,1)=2.0366806792244314E-003
uqx(3,2)=-3.0623567467657864E-003
uqy(3,2)=0.0000000000000000
uqx(4,1)=-4.4036485598978421E-004
uqy(4,1)=2.5646772695932669E-004
uqx(4,2)=-4.8184834784040605E-004
uqy(4,2)=1.1970911652618624E-004
uqx(4,3)=-4.8184834784040605E-004
uqy(4,3)=-1.1970911652618624E-004
uqx(5,1)=-6.4043997656082756E-005
uqy(5,1)=3.7294216529563565E-005
uqx(5,2)=-8.4819889884739453E-005
uqy(5,2)=2.9733239460941308E-005
uqx(5,3)=-8.2733545432358833E-005
uqy(5,3)=0.0000000000000000
uqx(5,4)=-8.4819889884739453E-005
uqy(5,4)=-2.9733239460941308E-005
endif
if(ntheta.EQ.31)then
uqx(2,1)=-0.04503436058954519
uqy(2,1)=0.026199536782378557
uqx(3,1)=-0.01026859419341395
uqy(3,1)=0.005985591626852251
uqx(3,2)=-0.005167330972412901
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0022619623413197005
uqy(4,1)=0.001318461208323861
uqx(4,2)=-0.0016192961158181562
uqy(4,2)=0.0004555838566847598
uqx(4,3)=-0.0016192961158181562
uqy(4,3)=-0.0004555838566847598
uqx(5,1)=-0.0005624625422492005
uqy(5,1)=0.00032774751316138077
uqx(5,2)=-0.0005173294370693253
uqy(5,2)=0.00019265284061534483
uqx(5,3)=-0.0004602795050385721
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0005173294370693253
uqy(5,4)=-0.00019265284061534483
endif
if(ntheta.EQ.32)then
uqx(2,1)=-0.01654760141946644
uqy(2,1)=0.00966452567895123
uqx(3,1)=-0.012359871245682334
uqy(3,1)=0.007207714818054002
uqx(3,2)=-0.003911727256791495
uqy(3,2)=0.0_dp
uqx(4,1)=-0.002858786939931005
uqy(4,1)=0.0016674960873579786
uqx(4,2)=-0.0016882532090136662
uqy(4,2)=0.0006763448089550626
uqx(4,3)=-0.0016882532090136662
uqy(4,3)=-0.0006763448089550626
uqx(5,1)=-0.0007371897937859137
uqy(5,1)=0.0004297244763339039
uqx(5,2)=-0.000627631582573753
uqy(5,2)=0.0002643079420545157
uqx(5,3)=-0.00048422193188717426
uqy(5,3)=0.0_dp
uqx(5,4)=-0.000627631582573753
uqy(5,4)=-0.0002643079420545157
endif
if(ntheta.EQ.33)then
uqx(2,1)=0.022424261931336693
uqy(2,1)=-0.01293583010556047
uqx(3,1)=-0.015051111904512509
uqy(3,1)=0.008780868998581349
uqx(3,2)=-0.0017190158015901534
uqy(3,2)=0.0_dp
uqx(4,1)=-0.003691155602273118
uqy(4,1)=0.002154674287298329
uqx(4,2)=-0.0014452371144368403
uqy(4,2)=0.0010001341858090018
uqx(4,3)=-0.0014452371144368403
uqy(4,3)=-0.0010001341858090018
uqx(5,1)=-0.0010068748692384517
uqy(5,1)=0.0005870299125574267
uqx(5,2)=-0.0007653530744728214
uqy(5,2)=0.00036152658483681714
uqx(5,3)=-0.0004553499125495603
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0007653530744728214
uqy(5,4)=-0.00036152658483681714
endif
if(ntheta.EQ.34)then
uqx(2,1)=0.04910578041208702
uqy(2,1)=-0.028397936546263408
uqx(3,1)=-0.017178989759262957
uqy(3,1)=0.01002246554565527
uqx(3,2)=8.264575740824117e-05
uqy(3,2)=0.0_dp
uqx(4,1)=-0.004400679112276507
uqy(4,1)=0.0025693001210308855
uqx(4,2)=-0.0009881999167575153
uqy(4,2)=0.0012622803233798047
uqx(4,3)=-0.0009881999167575153
uqy(4,3)=-0.0012622803233798047
uqx(5,1)=-0.001272219388203348
uqy(5,1)=0.0007415449672285427
uqx(5,2)=-0.000879355805314842
uqy(5,2)=0.0004375414558140443
uqx(5,3)=-0.00042505551681612546
uqy(5,3)=0.0_dp
uqx(5,4)=-0.000879355805314842
uqy(5,4)=-0.0004375414558140443
endif
if(ntheta.EQ.35)then
uqx(2,1)=0.06865399041348592
uqy(2,1)=-0.03972106739590168
uqx(3,1)=-0.018952730759702035
uqy(3,1)=0.011055322836226759
uqx(3,2)=0.0016473822918342902
uqy(3,2)=0.0_dp
uqx(4,1)=-0.005036952560980132
uqy(4,1)=0.0029403617530366715
uqx(4,2)=-0.00042017634429974027
uqy(4,2)=0.0014901107428717552
uqx(4,3)=-0.00042017634429974027
uqy(4,3)=-0.0014901107428717552
uqx(5,1)=-0.0015365634357938243
uqy(5,1)=0.0008952889713287348
uqx(5,2)=-0.0009801166864731095
uqy(5,2)=0.0005010813871941349
uqx(5,3)=-0.00040832163560121045
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0009801166864731095
uqy(5,4)=-0.0005010813871941349
endif

!!!!!!!!!!!!!!!!!!! Layer 1 at z=0

 do i=1,ndim/2

   rx=coord(i,1)
   ry=coord(i,2)

   delta(1)=(0.d0,0.d0)
   delta(2)=(0.d0,0.d0)

!write(*,*) n,m,vk(1),vk(2),rx,ry

   ncount=0
     do nn=1,5,2
       do n=2,5
         do m=1,n-1

           vk(1)=(n-1)*Gn(nn,1)+(m-1)*Gn(nn+2,1)
           vk(2)=(n-1)*Gn(nn,2)+(m-1)*Gn(nn+2,2)

           ux=cos((nn-1)*pi/3.0_dp)*uqx(n,m)-sin((nn-1)*pi/3.0_dp)*uqy(n,m)
           uy=sin((nn-1)*pi/3.0_dp)*uqx(n,m)+cos((nn-1)*pi/3.0_dp)*uqy(n,m)

           delta(1)=delta(1)-2*ux*sin(vk(1)*rx+vk(2)*ry)
           delta(2)=delta(2)-2*uy*sin(vk(1)*rx+vk(2)*ry)

         end do
       end do
     end do

!     if(coord(i,3).LT.deltaR) then
       coord(i,1)=rx+delta(1)/2
       coord(i,2)=ry+delta(2)/2
!     else


    coord(i+ndim/2,1)=-coord(i,1)!rx-delta(1)/2
    coord(i+ndim/2,2)=coord(i,2) !ry-delta(2)/2


enddo


end subroutine Relaxation_old

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Relaxation(coord,Gn)

 real(dp), intent(inout) :: coord(ndim,3)
 real(dp), intent(in)    :: Gn(8,2)

 integer(dp) :: i,j, m, n, ncount, nn
 real(dp) :: coordnr(ndim,3), coordtemp(ndim,3)
 real(dp) :: ux, uy, rx, ry
 real(dp) :: uqx(5,5), uqy(5,5), delta(2)
 real(dp) :: vk(2)
 real(dp), allocatable :: uq(:,:)

 allocate(uq(ndim,2))

 coordnr(:,:) = coord(:,:)

  uqx=0.0_dp
  uqy=0.0_dp

  if(ntheta.EQ.2)then
  uqx(2,1)=-0.0007312154977677281
  uqy(2,1)=0.0004714203951774816
  uqx(3,1)=-1.154897989107238e-06
  uqy(3,1)=7.447553144970968e-07
  uqx(3,2)=-1.586938089071483e-06
  uqy(3,2)=0.0_dp
  uqx(4,1)=-2.434676944887857e-09
  uqy(4,1)=1.570103901518326e-09
  uqx(4,2)=-4.341715993276229e-09
  uqy(4,2)=1.000979455471334e-09
  uqx(4,3)=-4.341715993276229e-09
  uqy(4,3)=-1.000979455471334e-09
  uqx(5,1)=-5.701274058766884e-12
  uqy(5,1)=3.740789872160164e-12
  uqx(5,2)=-1.288484816335295e-11
  uqy(5,2)=4.352172282562038e-12
  uqx(5,3)=-1.480922507260647e-11
  uqy(5,3)=0.0_dp
  uqx(5,4)=-1.288484816335295e-11
  uqy(5,4)=-4.352172282562038e-12
  endif
  if(ntheta.EQ.3)then
  uqx(2,1)=-0.001435896542049262
  uqy(2,1)=0.0008971722346957707
  uqx(3,1)=-4.431771204029964e-06
  uqy(3,1)=2.770018809111695e-06
  uqx(3,2)=-5.980337522539451e-06
  uqy(3,2)=0.0_dp
  uqx(4,1)=-1.828082429017813e-08
  uqy(4,1)=1.142651350640817e-08
  uqx(4,2)=-3.192212955644209e-08
  uqy(4,2)=7.41532195221355e-09
  uqx(4,3)=-3.192212955644209e-08
  uqy(4,3)=-7.41532195221355e-09
  uqx(5,1)=-8.489859929513374e-11
  uqy(5,1)=5.318288726255601e-11
  uqx(5,2)=-1.864781523972493e-10
  uqy(5,2)=6.321663026418908e-11
  uqx(5,3)=-2.126102770742755e-10
  uqy(5,3)=0.0_dp
  uqx(5,4)=-1.864781523972493e-10
  uqy(5,4)=-6.321663026418908e-11
  endif
  if(ntheta.EQ.4)then
  uqx(2,1)=-0.002370539233839467
  uqy(2,1)=0.001455154743491223
  uqx(3,1)=-1.207488833535449e-05
  uqy(3,1)=7.415571455420628e-06
  uqx(3,2)=-1.606356439211949e-05
  uqy(3,2)=0.0_dp
  uqx(4,1)=-8.232092034607089e-08
  uqy(4,1)=5.05576290960723e-08
  uqx(4,2)=-1.412786099990225e-07
  uqy(4,2)=3.296261058695678e-08
  uqx(4,3)=-1.412786099990225e-07
  uqy(4,3)=-3.296261058695678e-08
  uqx(5,1)=-6.338434146746167e-10
  uqy(5,1)=3.894401693009466e-10
  uqx(5,2)=-1.363506478716756e-09
  uqy(5,2)=4.635827193166175e-10
  uqx(5,3)=-1.545001808355566e-09
  uqy(5,3)=0.0_dp
  uqx(5,4)=-1.363506478716756e-09
  uqy(5,4)=-4.635827193166175e-10
  endif
  if(ntheta.EQ.5)then
  uqx(2,1)=-0.00352744756912305
  uqy(2,1)=0.002140730764691146
  uqx(3,1)=-2.68010445377339e-05
  uqy(3,1)=1.627423402644026e-05
  uqx(3,2)=-3.516363993145865e-05
  uqy(3,2)=0.0_dp
  uqx(4,1)=-2.729706294022432e-07
  uqy(4,1)=1.657616285518053e-07
  uqx(4,2)=-4.604302998995516e-07
  uqy(4,2)=1.078115792758269e-07
  uqx(4,3)=-4.604302998995516e-07
  uqy(4,3)=-1.078115792758269e-07
  uqx(5,1)=-3.147029295689003e-09
  uqy(5,1)=1.911153328301163e-09
  uqx(5,2)=-6.629538953061466e-09
  uqy(5,2)=2.259053162518054e-09
  uqx(5,3)=-7.46551582954267e-09
  uqy(5,3)=0.0_dp
  uqx(5,4)=-6.629538953061466e-09
  uqy(5,4)=-2.259053162518054e-09
  endif
  if(ntheta.EQ.6)then
  uqx(2,1)=-0.004895854840040079
  uqy(2,1)=0.002947504859599212
  uqx(3,1)=-5.185833778140336e-05
  uqy(3,1)=3.124215464415819e-05
  uqx(3,2)=-6.704921137114243e-05
  uqy(3,2)=0.0_dp
  uqx(4,1)=-7.375225497266297e-07
  uqy(4,1)=4.443469426422899e-07
  uqx(4,2)=-1.22125720193429e-06
  uqy(4,2)=2.869558314022102e-07
  uqx(4,3)=-1.22125720193429e-06
  uqy(4,3)=-2.869558314022102e-07
  uqx(5,1)=-1.189937627779894e-08
  uqy(5,1)=7.168907495586338e-09
  uqx(5,2)=-2.451139179182853e-08
  uqy(5,2)=8.36909403937007e-09
  uqx(5,3)=-2.741742519309818e-08
  uqy(5,3)=0.0_dp
  uqx(5,4)=-2.451139179182853e-08
  uqy(5,4)=-8.36909403937007e-09
  endif
  if(ntheta.EQ.7)then
  uqx(2,1)=-0.006461449348026135
  uqy(2,1)=0.00386704590013867
  uqx(3,1)=-9.088953324281387e-05
  uqy(3,1)=5.443895214647056e-05
  uqx(3,2)=-0.000115650920643144
  uqy(3,2)=0.0_dp
  uqx(4,1)=-1.719155170137875e-06
  uqy(4,1)=1.029772845664776e-06
  uqx(4,2)=-2.790313685737744e-06
  uqy(4,2)=6.580010547012858e-07
  uqx(4,3)=-2.790313685737744e-06
  uqy(4,3)=-6.580010547012858e-07
  uqx(5,1)=-3.697520081119331e-08
  uqy(5,1)=2.214633072942774e-08
  uqx(5,2)=-7.434248153831823e-08
  uqy(5,2)=2.543295223277813e-08
  uqx(5,3)=-8.254619477080724e-08
  uqy(5,3)=0.0_dp
  uqx(5,4)=-7.434248153831823e-08
  uqy(5,4)=-2.543295223277813e-08
  endif
  if(ntheta.EQ.8)then
  uqx(2,1)=-0.008206072919252466
  uqy(2,1)=0.004888718057559328
  uqx(3,1)=-0.0001477464619472436
  uqy(3,1)=8.80996557824374e-05
  uqx(3,2)=-0.000184730257310846
  uqy(3,2)=0.0_dp
  uqx(4,1)=-3.581496930090418e-06
  uqy(4,1)=2.135790438804529e-06
  uqx(4,2)=-5.688326741335424e-06
  uqy(4,2)=1.34655485744528e-06
  uqx(4,3)=-5.688326741335424e-06
  uqy(4,3)=-1.34655485744528e-06
  uqx(5,1)=-9.895113648508195e-08
  uqy(5,1)=5.900265937736671e-08
  uqx(5,2)=-1.938529162557286e-07
  uqy(5,2)=6.645103335141035e-08
  uqx(5,3)=-2.135258580530001e-07
  uqy(5,3)=0.0_dp
  uqx(5,4)=-1.938529162557286e-07
  uqy(5,4)=-6.645103335141035e-08
  endif
  if(ntheta.EQ.9)then
  uqx(2,1)=-0.01010772196045918
  uqy(2,1)=0.005999689628423322
  uqx(3,1)=-0.0002262612975656868
  uqy(3,1)=0.0001344419648152459
  uqx(3,2)=-0.0002775244610317592
  uqy(3,2)=0.0_dp
  uqx(4,1)=-6.825276745678818e-06
  uqy(4,1)=4.055911695604714e-06
  uqx(4,2)=-1.059052135909834e-05
  uqy(4,2)=2.517294978569754e-06
  uqx(4,3)=-1.059052135909834e-05
  uqy(4,3)=-2.517294978569754e-06
  uqx(5,1)=-2.352002772937919e-07
  uqy(5,1)=1.397520149575992e-07
  uqx(5,2)=-4.48274468435581e-07
  uqy(5,2)=1.539852461003042e-07
  uqx(5,3)=-4.895322080001663e-07
  uqy(5,3)=0.0_dp
  uqx(5,4)=-4.48274468435581e-07
  uqy(5,4)=-1.539852461003042e-07
  endif
  if(ntheta.EQ.10)then
  uqx(2,1)=-0.01214094723296179
  uqy(2,1)=0.007185172779869293
  uqx(3,1)=-0.0003299929115003222
  uqy(3,1)=0.000195519493609438
  uqx(3,2)=-0.0003964101199496858
  uqy(3,2)=0.0_dp
  uqx(4,1)=-1.209151294434386e-05
  uqy(4,1)=7.165004163245931e-06
  uqx(4,2)=-1.830230180218659e-05
  uqy(4,2)=4.369241246891779e-06
  uqx(4,3)=-1.830230180218659e-05
  uqy(4,3)=-4.369241246891779e-06
  uqx(5,1)=-5.072472977662023e-07
  uqy(5,1)=3.005414175984017e-07
  uqx(5,2)=-9.393528199777028e-07
  uqy(5,2)=3.233772873927477e-07
  uqx(5,3)=-1.016494567114082e-06
  uqy(5,3)=0.0_dp
  uqx(5,4)=-9.393528199777028e-07
  uqy(5,4)=-3.233772873927477e-07
  endif
  if(ntheta.EQ.11)then
  uqx(2,1)=-0.01427768000921253
  uqy(2,1)=0.008428908168456773
  uqx(3,1)=-0.0004619778319639407
  uqy(3,1)=0.0002730781672413903
  uqx(3,2)=-0.0005426360143594519
  uqy(3,2)=0.0_dp
  uqx(4,1)=-2.014728185682702e-05
  uqy(4,1)=1.191074841386241e-05
  uqx(4,2)=-2.970961882633886e-05
  uqy(4,2)=7.124796626544694e-06
  uqx(4,3)=-2.970961882633886e-05
  uqy(4,3)=-7.124796626544694e-06
  uqx(5,1)=-1.008004944770031e-06
  uqy(5,1)=5.958424731241988e-07
  uqx(5,2)=-1.811980725955097e-06
  uqy(5,2)=6.251978084599946e-07
  uqx(5,3)=-1.942211173240885e-06
  uqy(5,3)=0.0_dp
  uqx(5,4)=-1.811980725955097e-06
  uqy(5,4)=-6.251978084599946e-07
  endif
  if(ntheta.EQ.12)then
  uqx(2,1)=-0.01648841976833508
  uqy(2,1)=0.00971385525211049
  uqx(3,1)=-0.0006245208306885359
  uqy(3,1)=0.0003684360208645578
  uqx(3,2)=-0.0007161707066312485
  uqy(3,2)=0.0_dp
  uqx(4,1)=-3.185339267975134e-05
  uqy(4,1)=1.879461360555224e-05
  uqx(4,2)=-4.570953565691623e-05
  uqy(4,2)=1.101356078295402e-05
  uqx(4,3)=-4.570953565691623e-05
  uqy(4,3)=-1.101356078295402e-05
  uqx(5,1)=-1.86730786016839e-06
  uqy(5,1)=1.101634412158447e-06
  uqx(5,2)=-3.256105072314334e-06
  uqy(5,2)=1.12609910984101e-06
  uqx(5,3)=-3.456108802116457e-06
  uqy(5,3)=0.0_dp
  uqx(5,4)=-3.256105072314334e-06
  uqy(5,4)=-1.12609910984101e-06
  endif
  if(ntheta.EQ.13)then
  uqx(2,1)=-0.01874363020388613
  uqy(2,1)=0.01102299870360727
  uqx(3,1)=-0.000819056300796783
  uqy(3,1)=0.0004824043177345613
  uqx(3,2)=-0.0009156910682583513
  uqy(3,2)=0.0_dp
  uqx(4,1)=-4.81174507641687e-05
  uqy(4,1)=2.834441813303126e-05
  uqx(4,2)=-6.713199189077586e-05
  uqy(4,2)=1.62532984542057e-05
  uqx(4,3)=-6.713199189077586e-05
  uqy(4,3)=-1.62532984542057e-05
  uqx(5,1)=-3.254182136977055e-06
  uqy(5,1)=1.916674289995359e-06
  uqx(5,2)=-5.502334381872796e-06
  uqy(5,2)=1.907490988568425e-06
  uqx(5,3)=-5.782489009202849e-06
  uqy(5,3)=0.0_dp
  uqx(5,4)=-5.502334381872796e-06
  uqy(5,4)=-1.907490988568425e-06
  endif
  if(ntheta.EQ.14)then
  uqx(2,1)=-0.02101513878443223
  uqy(2,1)=0.01234015250707122
  uqx(3,1)=-0.001046098478254697
  uqy(3,1)=0.0006152602359553973
  uqx(3,2)=-0.001138709526361856
  uqy(3,2)=0.0_dp
  uqx(4,1)=-6.983936154117285e-05
  uqy(4,1)=4.108260075599009e-05
  uqx(4,2)=-9.46663508367406e-05
  uqy(4,2)=2.303123144864916e-05
  uqx(4,3)=-9.46663508367406e-05
  uqy(4,3)=-2.303123144864916e-05
  uqx(5,1)=-5.374846134548899e-06
  uqy(5,1)=3.161278639182164e-06
  uqx(5,2)=-8.811118519862508e-06
  uqy(5,2)=3.06194712722794e-06
  uqx(5,3)=-9.167731209069543e-06
  uqy(5,3)=0.0_dp
  uqx(5,4)=-8.811118519862508e-06
  uqy(5,4)=-3.06194712722794e-06
  endif
  if(ntheta.EQ.15)then
  uqx(2,1)=-0.0232773402036116
  uqy(2,1)=0.01365064721466022
  uqx(3,1)=-0.001305280475347721
  uqy(3,1)=0.0007667708872687428
  uqx(3,2)=-0.001381810113794881
  uqy(3,2)=0.0_dp
  uqx(4,1)=-9.785794354653753e-05
  uqy(4,1)=5.749524297849592e-05
  uqx(4,2)=-0.0001288049100968271
  uqy(4,2)=3.148869813455972e-05
  uqx(4,3)=-0.0001288049100968271
  uqy(4,3)=-3.148869813455972e-05
  uqx(5,1)=-8.46637610884184e-06
  uqy(5,1)=4.973588375291042e-06
  uqx(5,2)=-1.345696760939882e-05
  uqy(5,2)=4.687804727529729e-06
  uqx(5,3)=-1.386349319225168e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-1.345696760939882e-05
  uqy(5,4)=-4.687804727529729e-06
  endif
  if(ntheta.EQ.16)then
  uqx(2,1)=-0.02550806137259484
  uqy(2,1)=0.01494181959537789
  uqx(3,1)=-0.001595465429771206
  uqy(3,1)=0.000936258825655531
  uqx(3,2)=-0.001640945873196756
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0001329064560728502
  uqy(4,1)=7.80063592484412e-05
  uqx(4,2)=-0.0001698108360147586
  uqy(4,2)=4.171126842036178e-05
  uqx(4,3)=-0.0001698108360147586
  uqy(4,3)=-4.171126842036178e-05
  uqx(5,1)=-1.278694321771364e-05
  uqy(5,1)=7.503853321298472e-06
  uqx(5,2)=-1.971030956794678e-05
  uqy(5,2)=6.882830706202828e-06
  uqx(5,3)=-2.010885042278717e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-1.971030956794678e-05
  uqy(5,4)=-6.882830706202828e-06
  endif
  if(ntheta.EQ.17)then
  uqx(2,1)=-0.02768902996623525
  uqy(2,1)=0.01620327280964095
  uqx(3,1)=-0.001914903012504294
  uqy(3,1)=0.001122693458825686
  uqx(3,2)=-0.001911747148401908
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0001755820633807187
  uqy(4,1)=0.0001029603397676965
  uqx(4,2)=-0.0002177119608478438
  uqy(4,2)=5.372503557832448e-05
  uqx(4,3)=-0.0002177119608478438
  uqy(4,3)=-5.372503557832448e-05
  uqx(5,1)=-1.860418596133933e-05
  uqy(5,1)=1.090765125005818e-05
  uqx(5,2)=-2.78198841681919e-05
  uqy(5,2)=9.737957708444016e-06
  uqx(5,3)=-2.811427798376512e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-2.78198841681919e-05
  uqy(5,4)=-9.737957708444016e-06
  endif
  if(ntheta.EQ.18)then
  uqx(2,1)=-0.0298059688900068
  uqy(2,1)=0.01742692140290434
  uqx(3,1)=-0.002261402797319784
  uqy(3,1)=0.001324791852622696
  uqx(3,2)=-0.00218979975575966
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0002263307798910572
  uqy(4,1)=0.0001326134069799172
  uqx(4,2)=-0.0002723167306780222
  uqy(4,2)=6.749853746454043e-05
  uqx(4,3)=-0.0002723167306780222
  uqy(4,3)=-6.749853746454043e-05
  uqx(5,1)=-2.618344775700557e-05
  uqy(5,1)=1.533904799393413e-05
  uqx(5,2)=-3.799806596259685e-05
  uqy(5,2)=1.333194246897647e-05
  uqx(5,3)=-3.804957109168239e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-3.799806596259685e-05
  uqy(5,4)=-1.333194246897647e-05
  endif
  if(ntheta.EQ.19)then
  uqx(2,1)=-0.03184839354005421
  uqy(2,1)=0.01860686639767005
  uqx(3,1)=-0.002632500871417905
  uqy(3,1)=0.001541115305242353
  uqx(3,2)=-0.00247086793286604
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0002854464623341651
  uqy(4,1)=0.0001671332222719948
  uqx(4,2)=-0.0003332456229776692
  uqy(4,2)=8.294893631494263e-05
  uqx(4,3)=-0.0003332456229776692
  uqy(4,3)=-8.294893631494263e-05
  uqx(5,1)=-3.5777331916109e-05
  uqy(5,1)=2.094453998211611e-05
  uqx(5,2)=-5.041053364039542e-05
  uqy(5,2)=1.772747608189362e-05
  uqx(5,3)=-5.003665347329104e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-5.041053364039542e-05
  uqy(5,4)=-1.772747608189362e-05
  endif
  if(ntheta.EQ.20)then
  uqx(2,1)=-0.03380921098354144
  uqy(2,1)=0.01973915809034029
  uqx(3,1)=-0.003025604178510975
  uqy(3,1)=0.001770152796829164
  uqx(3,2)=-0.002751053010400932
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0003530806004474008
  uqy(4,1)=0.0002066047359070089
  uqx(4,2)=-0.0003999707148415089
  uqy(4,2)=9.995079644234232e-05
  uqx(4,3)=-0.0003999707148415089
  uqy(4,3)=-9.995079644234232e-05
  uqx(5,1)=-4.761749968378625e-05
  uqy(5,1)=2.785831220463549e-05
  uqx(5,2)=-6.517066822235473e-05
  uqy(5,2)=2.296892009302885e-05
  uqx(5,3)=-6.414716777612842e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-6.517066822235473e-05
  uqy(5,4)=-2.296892009302885e-05
  endif
  if(ntheta.EQ.21)then
  uqx(2,1)=-0.03568421614472853
  uqy(2,1)=0.0208215014107083
  uqx(3,1)=-0.003438105278930305
  uqy(3,1)=0.002010387168015211
  uqx(3,2)=-0.003026890903307585
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0004292590682540321
  uqy(4,1)=0.0002510400416038367
  uqx(4,2)=-0.0004718571203498731
  uqy(4,2)=0.0001183459294007554
  uqx(4,3)=-0.0004718571203498731
  uqy(4,3)=-0.0001183459294007554
  uqx(5,1)=-6.190907493630538e-05
  uqy(5,1)=3.61990173774744e-05
  uqx(5,2)=-8.233825903866238e-05
  uqy(5,2)=2.908155778293243e-05
  uqx(5,3)=-8.040403158689358e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-8.233825903866238e-05
  uqy(5,4)=-2.908155778293243e-05
  endif
  if(ntheta.EQ.22)then
  uqx(2,1)=-0.03747156089826429
  uqy(2,1)=0.02185294737277631
  uqx(3,1)=-0.003867466589799954
  uqy(3,1)=0.002260343560684938
  uqx(3,2)=-0.003295398510176307
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0005139023086210218
  uqy(4,1)=0.0003003901861872707
  uqx(4,2)=-0.0005482018340347233
  uqy(4,2)=0.0001379531343067409
  uqx(4,3)=-0.0005482018340347233
  uqy(4,3)=-0.0001379531343067409
  uqx(5,1)=-7.882756813514927e-05
  uqy(5,1)=4.606802067573352e-05
  uqx(5,2)=-0.0001019216358938264
  uqy(5,2)=3.607207962412092e-05
  uqx(5,3)=-9.878582619642953e-05
  uqy(5,3)=0.0_dp
  uqx(5,4)=-0.0001019216358938264
  uqy(5,4)=-3.607207962412092e-05
  endif
  if(ntheta.EQ.23)then
  uqx(2,1)=-0.03917124800252712
  uqy(2,1)=0.02283360020288869
  uqx(3,1)=-0.00431127715331047
  uqy(3,1)=0.002518621941051702
  uqx(3,2)=-0.003554081631682294
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0006068462026936415
  uqy(4,1)=0.0003545573425025866
  uqx(4,2)=-0.0006282673934951395
  uqy(4,2)=0.00015857707553736
  uqx(4,3)=-0.0006282673934951395
  uqy(4,3)=-0.00015857707553736
  uqx(5,1)=-9.851794930567258e-05
  uqy(5,1)=5.754889134915976e-05
  uqx(5,2)=-0.0001238821974354894
  uqy(5,2)=4.392995908857313e-05
  uqx(5,3)=-0.0001192328797561733
  uqy(5,3)=0.0_dp
  uqx(5,4)=-0.0001238821974354894
  uqy(5,4)=-4.392995908857313e-05
  endif
  if(ntheta.EQ.24)then
  uqx(2,1)=-0.04078467983132041
  uqy(2,1)=0.02376435704408403
  uqx(3,1)=-0.004767286855073334
  uqy(3,1)=0.00278391658984503
  uqx(3,2)=-0.003800916556729175
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0007078617681184562
  uqy(4,1)=0.0004134062745389956
  uqx(4,2)=-0.0007113093010572949
  uqy(4,2)=0.0001800159036691447
  uqx(4,3)=-0.0007113093010572949
  uqy(4,3)=-0.0001800159036691447
  uqx(5,1)=-0.0001210953749828021
  uqy(5,1)=7.070785140253987e-05
  uqx(5,2)=-0.0001481403717251601
  uqy(5,2)=5.26293883724199e-05
  uqx(5,3)=-0.0001416540841691958
  uqy(5,3)=0.0_dp
  uqx(5,4)=-0.0001481403717251601
  uqy(5,4)=-5.26293883724199e-05
  endif
  if(ntheta.EQ.25)then
  uqx(2,1)=-0.04231427495534361
  uqy(2,1)=0.02464668743107288
  uqx(3,1)=-0.005233423406853966
  uqy(3,1)=0.003055025655109403
  uqx(3,2)=-0.004034315476632931
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.000816672621751051
  uqy(4,1)=0.0004767744840141218
  uqx(4,2)=-0.0007965971920023176
  uqy(4,2)=0.0002020674940364024
  uqx(4,3)=-0.0007965971920023176
  uqy(4,3)=-0.0002020674940364024
  uqx(5,1)=-0.0001466470646981058
  uqy(5,1)=8.559488727257606e-05
  uqx(5,2)=-0.0001745822216037299
  uqy(5,2)=6.213149817457247e-05
  uqx(5,3)=-0.0001659337332567695
  uqy(5,3)=0.0_dp
  uqx(5,4)=-0.0001745822216037299
  uqy(5,4)=-6.213149817457247e-05
  endif
  if(ntheta.EQ.26)then
  uqx(2,1)=-0.0437631542688065
  uqy(2,1)=0.02548245325111225
  uqx(3,1)=-0.005707796948781806
  uqy(3,1)=0.003330853585610454
  uqx(3,2)=-0.00425308343192316
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.0009329697453595618
  uqy(4,1)=0.000544480777040355
  uqx(4,2)=-0.0008834303440076407
  uqy(4,2)=0.0002245343486965423
  uqx(4,3)=-0.0008834303440076407
  uqy(4,3)=-0.0002245343486965423
  uqx(5,1)=-0.0001752348810775276
  uqy(5,1)=0.0001022452652849156
  uqx(5,2)=-0.0002030661159222708
  uqy(5,2)=7.23866540137963e-05
  uqx(5,3)=-0.0001919379139721673
  uqy(5,3)=0.0_dp
  uqx(5,4)=-0.0002030661159222708
  uqy(5,4)=-7.23866540137963e-05
  endif
  if(ntheta.EQ.27)then
  uqx(2,1)=-0.04513489169976175
  uqy(2,1)=0.02627376613853282
  uqx(3,1)=-0.006188696292743292
  uqy(3,1)=0.003610408772998517
  uqx(3,2)=-0.0044563721863212
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.001056423499571107
  uqy(4,1)=0.0006163322244658482
  uqx(4,2)=-0.0009711483975435052
  uqy(4,2)=0.0002472273003311512
  uqx(4,3)=-0.0009711483975435052
  uqy(4,3)=-0.0002472273003311512
  uqx(5,1)=-0.0002068982555902986
  uqy(5,1)=0.0001206812433280477
  uqx(5,2)=-0.0002334290810502318
  uqy(5,2)=8.333668679649319e-05
  uqx(5,3)=-0.0002195201853748941
  uqy(5,3)=0.0_dp
  uqx(5,4)=-0.0002334290810502318
  uqy(5,4)=-8.333668679649319e-05
  endif
  if(ntheta.EQ.28)then
  uqx(2,1)=-0.04643332085491404
  uqy(2,1)=0.02702287718328207
  uqx(3,1)=-0.006674579953450168
  uqy(3,1)=0.003892798220305804
  uqx(3,2)=-0.004643634524000381
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.001186693080114393
  uqy(4,1)=0.0006921296326732865
  uqx(4,2)=-0.001059138218690306
  uqy(4,2)=0.0002699681973246943
  uqx(4,3)=-0.001059138218690306
  uqy(4,3)=-0.0002699681973246943
  uqx(5,1)=-0.0002416571953113832
  uqy(5,1)=0.000140913825684749
  uqx(5,2)=-0.0002654926052230902
  uqy(5,2)=9.491697014352015e-05
  uqx(5,3)=-0.000248526431033852
  uqy(5,3)=0.0_dp
  uqx(5,4)=-0.0002654926052230902
  uqy(5,4)=-9.491697014352015e-05
  endif
  if(ntheta.EQ.29)then
  uqx(2,1)=-0.0476623722906518
  uqy(2,1)=0.02773208419857225
  uqx(3,1)=-0.007164065345384643
  uqy(3,1)=0.004177221193478502
  uqx(3,2)=-0.004814580805092504
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.001323433931627007
  uqy(4,1)=0.0007716718280408433
  uqx(4,2)=-0.001146837865060967
  uqy(4,2)=0.0002925918557552575
  uqx(4,3)=-0.001146837865060967
  uqy(4,3)=-0.0002925918557552575
  uqx(5,1)=-0.0002795152299791238
  uqy(5,1)=0.0001629444802374371
  uqx(5,2)=-0.0002990678227357772
  uqy(5,2)=0.0001070583234817996
  uqx(5,3)=-0.0002787988860265953
  uqy(5,3)=0.0_dp
  uqx(5,4)=-0.0002990678227357772
  uqy(5,4)=-0.0001070583234817996
  endif
  if(ntheta.EQ.30)then
  uqx(2,1)=-0.04882563370160641
  uqy(2,1)=0.02840347765698087
  uqx(3,1)=-0.00765593696089868
  uqy(3,1)=0.004462974022425572
  uqx(3,2)=-0.004969132784832358
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.00146630725510065
  uqy(4,1)=0.0008547611790123954
  uqx(4,2)=-0.001233739988517768
  uqy(4,2)=0.0003149494088037902
  uqx(4,3)=-0.001233739988517768
  uqy(4,3)=-0.0003149494088037902
  uqx(5,1)=-0.0003204630115092076
  uqy(5,1)=0.0001867672361436166
  uqx(5,2)=-0.0003339607774787402
  uqy(5,2)=0.0001196892313188401
  uqx(5,3)=-0.0003101796338676208
  uqy(5,3)=0.0_dp
  uqx(5,4)=-0.0003339607774787402
  uqy(5,4)=-0.0001196892313188401
  endif
  if(ntheta.EQ.31)then
  uqx(2,1)=-0.04992193062183704
  uqy(2,1)=0.0290363736961492
  uqx(3,1)=-0.008149400568775942
  uqy(3,1)=0.004749598900816575
  uqx(3,2)=-0.005107282626575934
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.001615038562818766
  uqy(4,1)=0.0009412378657909575
  uqx(4,2)=-0.001319412580829482
  uqy(4,2)=0.0003369375238227954
  uqx(4,3)=-0.001319412580829482
  uqy(4,3)=-0.0003369375238227954
  uqx(5,1)=-0.0003644932065908745
  uqy(5,1)=0.0002123773891528697
  uqx(5,2)=-0.0003699867939099062
  uqy(5,2)=0.0001327448853399043
  uqx(5,3)=-0.0003425164020803855
  uqy(5,3)=0.0_dp
  uqx(5,4)=-0.0003699867939099062
  uqy(5,4)=-0.0001327448853399043
  endif
  if(ntheta.EQ.32)then
  uqx(2,1)=-0.05090329799667861
  uqy(2,1)=0.02960289728720748
  uqx(3,1)=-0.008646542228826316
  uqy(3,1)=0.00503832273936899
  uqx(3,2)=-0.005227869327712846
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.001769947347072314
  uqy(4,1)=0.001031290132427384
  uqx(4,2)=-0.001403673062984919
  uqy(4,2)=0.0003587752109830483
  uqx(4,3)=-0.001403673062984919
  uqy(4,3)=-0.0003587752109830483
  uqx(5,1)=-0.0004117279745827106
  uqy(5,1)=0.0002398460686580247
  uqx(5,2)=-0.0004070748751161033
  uqy(5,2)=0.0001462439658497053
  uqx(5,3)=-0.0003756924781756205
  uqy(5,3)=0.0_dp
  uqx(5,4)=-0.0004070748751161033
  uqy(5,4)=-0.0001462439658497053
  endif
  if(ntheta.EQ.33)then
  uqx(2,1)=-0.0513606096830928
  uqy(2,1)=0.02986539511944403
  uqx(3,1)=-0.009170005042813718
  uqy(3,1)=0.005342398460657084
  uqx(3,2)=-0.005318370702859235
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.001935879097096427
  uqy(4,1)=0.001127757228004547
  uqx(4,2)=-0.001487689894477332
  uqy(4,2)=0.0003830788831587224
  uqx(4,3)=-0.001487689894477332
  uqy(4,3)=-0.0003830788831587224
  uqx(5,1)=-0.000463394046633049
  uqy(5,1)=0.0002698905192682427
  uqx(5,2)=-0.0004460132086709931
  uqy(5,2)=0.0001608687247518319
  uqx(5,3)=-0.0004098000837104587
  uqy(5,3)=0.0_dp
  uqx(5,4)=-0.0004460132086709931
  uqy(5,4)=-0.0001608687247518319
  endif
  if(ntheta.EQ.34)then
  uqx(2,1)=-0.04874144182744936
  uqy(2,1)=0.0283417826549473
  uqx(3,1)=-0.009858487467627138
  uqy(3,1)=0.005742969922777688
  uqx(3,2)=-0.005287936990165755
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.002144796524776807
  uqy(4,1)=0.001249360187141939
  uqx(4,2)=-0.001578285629772355
  uqy(4,2)=0.0004267732968133611
  uqx(4,3)=-0.001578285629772355
  uqy(4,3)=-0.0004267732968133611
  uqx(5,1)=-0.0005276562550504655
  uqy(5,1)=0.000307285033763345
  uqx(5,2)=-0.0004925301179564875
  uqy(5,2)=0.0001813270867116877
  uqx(5,3)=-0.0004456179148133969
  uqy(5,3)=0.0_dp
  uqx(5,4)=-0.0004925301179564875
  uqy(5,4)=-0.0001813270867116877
  endif
  if(ntheta.EQ.35)then
  uqx(2,1)=-0.03297296564955147
  uqy(2,1)=0.01918709821552975
  uqx(3,1)=-0.01123505415890843
  uqy(3,1)=0.006546108167901785
  uqx(3,2)=-0.004684706983380724
  uqy(3,2)=0.0_dp
  uqx(4,1)=-0.002533504540963353
  uqy(4,1)=0.001476247589033027
  uqx(4,2)=-0.001676066560316715
  uqy(4,2)=0.0005604350432240707
  uqx(4,3)=-0.001676066560316715
  uqy(4,3)=-0.0005604350432240707
  uqx(5,1)=-0.0006414180926895605
  uqy(5,1)=0.0003736003362540499
  uqx(5,2)=-0.0005690638663033982
  uqy(5,2)=0.0002272281056786763
  uqx(5,3)=-0.0004791697292710052
  uqy(5,3)=0.0_dp
  uqx(5,4)=-0.0005690638663033982
  uqy(5,4)=-0.0002272281056786763
  endif

  !!!!!!!!!!!!!!!!!!! Layer 1 at z=0

  do i=1,ndim/2

    rx=coord(i,1)
    ry=coord(i,2)

    delta(1)=(0.d0,0.d0)
    delta(2)=(0.d0,0.d0)

    !write(*,*) n,m,vk(1),vk(2),rx,ry
    do nn=1,5,2
      do n=2,5
        do m=1,n-1

          vk(1)=(n-1)*Gn(nn,1)+(m-1)*Gn(nn+2,1)
          vk(2)=(n-1)*Gn(nn,2)+(m-1)*Gn(nn+2,2)

          ux=cos((nn-1)*pi/3.0_dp)*uqx(n,m)-sin((nn-1)*pi/3.0_dp)*uqy(n,m)
          uy=sin((nn-1)*pi/3.0_dp)*uqx(n,m)+cos((nn-1)*pi/3.0_dp)*uqy(n,m)

          delta(1)=delta(1)-2*ux*sin(vk(1)*rx+vk(2)*ry)
          delta(2)=delta(2)-2*uy*sin(vk(1)*rx+vk(2)*ry)

        end do
      end do
    end do

    coord(i,1)=rx+delta(1)/2
    coord(i,2)=ry+delta(2)/2

    coord(i+ndim/2,1)=-coord(i,1)!rx-delta(1)/2
    coord(i+ndim/2,2)=coord(i,2) !ry-delta(2)/2

  enddo

  ncount=0
  do i=ndim/2+1,3*ndim/4
    do j=ndim/2+1,3*ndim/4
      if(sqrt((coord(i,1) - coordnr(j,1))**2 + (coord(i,2) - coordnr(j,2))**2).lt.0.25_dp)then
        ncount = ncount+1
        coordtemp(j,:) = coord(i,:)
        exit
      endif
    enddo
  enddo

  ncount=ncount+1
  coordtemp(3*ndim/4,1) = coord(3*ndim/4,1)
  coordtemp(3*ndim/4,2) = -coord(3*ndim/4,2)
  coordtemp(3*ndim/4,3) = coord(3*ndim/4,3)

  if(ncount.ne.ndim/4)then
    write(*,*) 'ERROR RELAXATION', ncount
  endif

  ncount=0
  do i=3*ndim/4+1,ndim
    do j=3*ndim/4+1,ndim
      if(sqrt((coord(i,1) - coordnr(j,1))**2 + (coord(i,2) - coordnr(j,2))**2).lt.0.25_dp)then
        ncount = ncount+1
        coordtemp(j,:) = coord(i,:)
        exit
      endif
    enddo
  enddo

  ncount=ncount+1
  coordtemp(ndim,1) = coord(ndim,1)
  coordtemp(ndim,2) = -coord(ndim,2)
  coordtemp(ndim,3) = coord(ndim,3)

  if(ncount.ne.ndim/4)then
    write(*,*) 'ERROR RELAXATION', ncount
  endif

  coord(ndim/2+1:ndim,:) = coordtemp(ndim/2+1:ndim,:)


end subroutine Relaxation



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Relaxation_2(coord,Gn)

 real(dp), intent(inout) :: coord(ndim,3)
 real(dp), intent(in)    :: Gn(8,2)

 integer(dp) :: i, m, n, ncount, nn
 real(dp) :: ux, uy, rx, ry
 real(dp) :: uqx(5,5), uqy(5,5), delta(2)
 real(dp) :: vk(2)
 real(dp), allocatable :: uq(:,:)

 allocate(uq(ndim,2))

uqx=0.0_dp
uqy=0.0_dp

if(ntheta.EQ.2)then
uqx(2,1)=-0.0007312154977677281
uqy(2,1)=0.0004714203951774816
uqx(3,1)=-1.154897989107238e-06
uqy(3,1)=7.447553144970968e-07
uqx(3,2)=-1.586938089071483e-06
uqy(3,2)=0.0_dp
uqx(4,1)=-2.434676944887857e-09
uqy(4,1)=1.570103901518326e-09
uqx(4,2)=-4.341715993276229e-09
uqy(4,2)=1.000979455471334e-09
uqx(4,3)=-4.341715993276229e-09
uqy(4,3)=-1.000979455471334e-09
uqx(5,1)=-5.701274058766884e-12
uqy(5,1)=3.740789872160164e-12
uqx(5,2)=-1.288484816335295e-11
uqy(5,2)=4.352172282562038e-12
uqx(5,3)=-1.480922507260647e-11
uqy(5,3)=0.0_dp
uqx(5,4)=-1.288484816335295e-11
uqy(5,4)=-4.352172282562038e-12
endif
if(ntheta.EQ.3)then
uqx(2,1)=-0.001435896542049262
uqy(2,1)=0.0008971722346957707
uqx(3,1)=-4.431771204029964e-06
uqy(3,1)=2.770018809111695e-06
uqx(3,2)=-5.980337522539451e-06
uqy(3,2)=0.0_dp
uqx(4,1)=-1.828082429017813e-08
uqy(4,1)=1.142651350640817e-08
uqx(4,2)=-3.192212955644209e-08
uqy(4,2)=7.41532195221355e-09
uqx(4,3)=-3.192212955644209e-08
uqy(4,3)=-7.41532195221355e-09
uqx(5,1)=-8.489859929513374e-11
uqy(5,1)=5.318288726255601e-11
uqx(5,2)=-1.864781523972493e-10
uqy(5,2)=6.321663026418908e-11
uqx(5,3)=-2.126102770742755e-10
uqy(5,3)=0.0_dp
uqx(5,4)=-1.864781523972493e-10
uqy(5,4)=-6.321663026418908e-11
endif
if(ntheta.EQ.4)then
uqx(2,1)=-0.002370539233839467
uqy(2,1)=0.001455154743491223
uqx(3,1)=-1.207488833535449e-05
uqy(3,1)=7.415571455420628e-06
uqx(3,2)=-1.606356439211949e-05
uqy(3,2)=0.0_dp
uqx(4,1)=-8.232092034607089e-08
uqy(4,1)=5.05576290960723e-08
uqx(4,2)=-1.412786099990225e-07
uqy(4,2)=3.296261058695678e-08
uqx(4,3)=-1.412786099990225e-07
uqy(4,3)=-3.296261058695678e-08
uqx(5,1)=-6.338434146746167e-10
uqy(5,1)=3.894401693009466e-10
uqx(5,2)=-1.363506478716756e-09
uqy(5,2)=4.635827193166175e-10
uqx(5,3)=-1.545001808355566e-09
uqy(5,3)=0.0_dp
uqx(5,4)=-1.363506478716756e-09
uqy(5,4)=-4.635827193166175e-10
endif
if(ntheta.EQ.5)then
uqx(2,1)=-0.00352744756912305
uqy(2,1)=0.002140730764691146
uqx(3,1)=-2.68010445377339e-05
uqy(3,1)=1.627423402644026e-05
uqx(3,2)=-3.516363993145865e-05
uqy(3,2)=0.0_dp
uqx(4,1)=-2.729706294022432e-07
uqy(4,1)=1.657616285518053e-07
uqx(4,2)=-4.604302998995516e-07
uqy(4,2)=1.078115792758269e-07
uqx(4,3)=-4.604302998995516e-07
uqy(4,3)=-1.078115792758269e-07
uqx(5,1)=-3.147029295689003e-09
uqy(5,1)=1.911153328301163e-09
uqx(5,2)=-6.629538953061466e-09
uqy(5,2)=2.259053162518054e-09
uqx(5,3)=-7.46551582954267e-09
uqy(5,3)=0.0_dp
uqx(5,4)=-6.629538953061466e-09
uqy(5,4)=-2.259053162518054e-09
endif
if(ntheta.EQ.6)then
uqx(2,1)=-0.004895854840040079
uqy(2,1)=0.002947504859599212
uqx(3,1)=-5.185833778140336e-05
uqy(3,1)=3.124215464415819e-05
uqx(3,2)=-6.704921137114243e-05
uqy(3,2)=0.0_dp
uqx(4,1)=-7.375225497266297e-07
uqy(4,1)=4.443469426422899e-07
uqx(4,2)=-1.22125720193429e-06
uqy(4,2)=2.869558314022102e-07
uqx(4,3)=-1.22125720193429e-06
uqy(4,3)=-2.869558314022102e-07
uqx(5,1)=-1.189937627779894e-08
uqy(5,1)=7.168907495586338e-09
uqx(5,2)=-2.451139179182853e-08
uqy(5,2)=8.36909403937007e-09
uqx(5,3)=-2.741742519309818e-08
uqy(5,3)=0.0_dp
uqx(5,4)=-2.451139179182853e-08
uqy(5,4)=-8.36909403937007e-09
endif
if(ntheta.EQ.7)then
uqx(2,1)=-0.006461449348026135
uqy(2,1)=0.00386704590013867
uqx(3,1)=-9.088953324281387e-05
uqy(3,1)=5.443895214647056e-05
uqx(3,2)=-0.000115650920643144
uqy(3,2)=0.0_dp
uqx(4,1)=-1.719155170137875e-06
uqy(4,1)=1.029772845664776e-06
uqx(4,2)=-2.790313685737744e-06
uqy(4,2)=6.580010547012858e-07
uqx(4,3)=-2.790313685737744e-06
uqy(4,3)=-6.580010547012858e-07
uqx(5,1)=-3.697520081119331e-08
uqy(5,1)=2.214633072942774e-08
uqx(5,2)=-7.434248153831823e-08
uqy(5,2)=2.543295223277813e-08
uqx(5,3)=-8.254619477080724e-08
uqy(5,3)=0.0_dp
uqx(5,4)=-7.434248153831823e-08
uqy(5,4)=-2.543295223277813e-08
endif
if(ntheta.EQ.8)then
uqx(2,1)=-0.008206072919252466
uqy(2,1)=0.004888718057559328
uqx(3,1)=-0.0001477464619472436
uqy(3,1)=8.80996557824374e-05
uqx(3,2)=-0.000184730257310846
uqy(3,2)=0.0_dp
uqx(4,1)=-3.581496930090418e-06
uqy(4,1)=2.135790438804529e-06
uqx(4,2)=-5.688326741335424e-06
uqy(4,2)=1.34655485744528e-06
uqx(4,3)=-5.688326741335424e-06
uqy(4,3)=-1.34655485744528e-06
uqx(5,1)=-9.895113648508195e-08
uqy(5,1)=5.900265937736671e-08
uqx(5,2)=-1.938529162557286e-07
uqy(5,2)=6.645103335141035e-08
uqx(5,3)=-2.135258580530001e-07
uqy(5,3)=0.0_dp
uqx(5,4)=-1.938529162557286e-07
uqy(5,4)=-6.645103335141035e-08
endif
if(ntheta.EQ.9)then
uqx(2,1)=-0.01010772196045918
uqy(2,1)=0.005999689628423322
uqx(3,1)=-0.0002262612975656868
uqy(3,1)=0.0001344419648152459
uqx(3,2)=-0.0002775244610317592
uqy(3,2)=0.0_dp
uqx(4,1)=-6.825276745678818e-06
uqy(4,1)=4.055911695604714e-06
uqx(4,2)=-1.059052135909834e-05
uqy(4,2)=2.517294978569754e-06
uqx(4,3)=-1.059052135909834e-05
uqy(4,3)=-2.517294978569754e-06
uqx(5,1)=-2.352002772937919e-07
uqy(5,1)=1.397520149575992e-07
uqx(5,2)=-4.48274468435581e-07
uqy(5,2)=1.539852461003042e-07
uqx(5,3)=-4.895322080001663e-07
uqy(5,3)=0.0_dp
uqx(5,4)=-4.48274468435581e-07
uqy(5,4)=-1.539852461003042e-07
endif
if(ntheta.EQ.10)then
uqx(2,1)=-0.01214094723296179
uqy(2,1)=0.007185172779869293
uqx(3,1)=-0.0003299929115003222
uqy(3,1)=0.000195519493609438
uqx(3,2)=-0.0003964101199496858
uqy(3,2)=0.0_dp
uqx(4,1)=-1.209151294434386e-05
uqy(4,1)=7.165004163245931e-06
uqx(4,2)=-1.830230180218659e-05
uqy(4,2)=4.369241246891779e-06
uqx(4,3)=-1.830230180218659e-05
uqy(4,3)=-4.369241246891779e-06
uqx(5,1)=-5.072472977662023e-07
uqy(5,1)=3.005414175984017e-07
uqx(5,2)=-9.393528199777028e-07
uqy(5,2)=3.233772873927477e-07
uqx(5,3)=-1.016494567114082e-06
uqy(5,3)=0.0_dp
uqx(5,4)=-9.393528199777028e-07
uqy(5,4)=-3.233772873927477e-07
endif
if(ntheta.EQ.11)then
uqx(2,1)=-0.01427768000921253
uqy(2,1)=0.008428908168456773
uqx(3,1)=-0.0004619778319639407
uqy(3,1)=0.0002730781672413903
uqx(3,2)=-0.0005426360143594519
uqy(3,2)=0.0_dp
uqx(4,1)=-2.014728185682702e-05
uqy(4,1)=1.191074841386241e-05
uqx(4,2)=-2.970961882633886e-05
uqy(4,2)=7.124796626544694e-06
uqx(4,3)=-2.970961882633886e-05
uqy(4,3)=-7.124796626544694e-06
uqx(5,1)=-1.008004944770031e-06
uqy(5,1)=5.958424731241988e-07
uqx(5,2)=-1.811980725955097e-06
uqy(5,2)=6.251978084599946e-07
uqx(5,3)=-1.942211173240885e-06
uqy(5,3)=0.0_dp
uqx(5,4)=-1.811980725955097e-06
uqy(5,4)=-6.251978084599946e-07
endif
if(ntheta.EQ.12)then
uqx(2,1)=-0.01648841976833508
uqy(2,1)=0.00971385525211049
uqx(3,1)=-0.0006245208306885359
uqy(3,1)=0.0003684360208645578
uqx(3,2)=-0.0007161707066312485
uqy(3,2)=0.0_dp
uqx(4,1)=-3.185339267975134e-05
uqy(4,1)=1.879461360555224e-05
uqx(4,2)=-4.570953565691623e-05
uqy(4,2)=1.101356078295402e-05
uqx(4,3)=-4.570953565691623e-05
uqy(4,3)=-1.101356078295402e-05
uqx(5,1)=-1.86730786016839e-06
uqy(5,1)=1.101634412158447e-06
uqx(5,2)=-3.256105072314334e-06
uqy(5,2)=1.12609910984101e-06
uqx(5,3)=-3.456108802116457e-06
uqy(5,3)=0.0_dp
uqx(5,4)=-3.256105072314334e-06
uqy(5,4)=-1.12609910984101e-06
endif
if(ntheta.EQ.13)then
uqx(2,1)=-0.01874363020388613
uqy(2,1)=0.01102299870360727
uqx(3,1)=-0.000819056300796783
uqy(3,1)=0.0004824043177345613
uqx(3,2)=-0.0009156910682583513
uqy(3,2)=0.0_dp
uqx(4,1)=-4.81174507641687e-05
uqy(4,1)=2.834441813303126e-05
uqx(4,2)=-6.713199189077586e-05
uqy(4,2)=1.62532984542057e-05
uqx(4,3)=-6.713199189077586e-05
uqy(4,3)=-1.62532984542057e-05
uqx(5,1)=-3.254182136977055e-06
uqy(5,1)=1.916674289995359e-06
uqx(5,2)=-5.502334381872796e-06
uqy(5,2)=1.907490988568425e-06
uqx(5,3)=-5.782489009202849e-06
uqy(5,3)=0.0_dp
uqx(5,4)=-5.502334381872796e-06
uqy(5,4)=-1.907490988568425e-06
endif
if(ntheta.EQ.14)then
uqx(2,1)=-0.02101513878443223
uqy(2,1)=0.01234015250707122
uqx(3,1)=-0.001046098478254697
uqy(3,1)=0.0006152602359553973
uqx(3,2)=-0.001138709526361856
uqy(3,2)=0.0_dp
uqx(4,1)=-6.983936154117285e-05
uqy(4,1)=4.108260075599009e-05
uqx(4,2)=-9.46663508367406e-05
uqy(4,2)=2.303123144864916e-05
uqx(4,3)=-9.46663508367406e-05
uqy(4,3)=-2.303123144864916e-05
uqx(5,1)=-5.374846134548899e-06
uqy(5,1)=3.161278639182164e-06
uqx(5,2)=-8.811118519862508e-06
uqy(5,2)=3.06194712722794e-06
uqx(5,3)=-9.167731209069543e-06
uqy(5,3)=0.0_dp
uqx(5,4)=-8.811118519862508e-06
uqy(5,4)=-3.06194712722794e-06
endif
if(ntheta.EQ.15)then
uqx(2,1)=-0.0232773402036116
uqy(2,1)=0.01365064721466022
uqx(3,1)=-0.001305280475347721
uqy(3,1)=0.0007667708872687428
uqx(3,2)=-0.001381810113794881
uqy(3,2)=0.0_dp
uqx(4,1)=-9.785794354653753e-05
uqy(4,1)=5.749524297849592e-05
uqx(4,2)=-0.0001288049100968271
uqy(4,2)=3.148869813455972e-05
uqx(4,3)=-0.0001288049100968271
uqy(4,3)=-3.148869813455972e-05
uqx(5,1)=-8.46637610884184e-06
uqy(5,1)=4.973588375291042e-06
uqx(5,2)=-1.345696760939882e-05
uqy(5,2)=4.687804727529729e-06
uqx(5,3)=-1.386349319225168e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-1.345696760939882e-05
uqy(5,4)=-4.687804727529729e-06
endif
if(ntheta.EQ.16)then
uqx(2,1)=-0.02550806137259484
uqy(2,1)=0.01494181959537789
uqx(3,1)=-0.001595465429771206
uqy(3,1)=0.000936258825655531
uqx(3,2)=-0.001640945873196756
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0001329064560728502
uqy(4,1)=7.80063592484412e-05
uqx(4,2)=-0.0001698108360147586
uqy(4,2)=4.171126842036178e-05
uqx(4,3)=-0.0001698108360147586
uqy(4,3)=-4.171126842036178e-05
uqx(5,1)=-1.278694321771364e-05
uqy(5,1)=7.503853321298472e-06
uqx(5,2)=-1.971030956794678e-05
uqy(5,2)=6.882830706202828e-06
uqx(5,3)=-2.010885042278717e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-1.971030956794678e-05
uqy(5,4)=-6.882830706202828e-06
endif
if(ntheta.EQ.17)then
uqx(2,1)=-0.02768902996623525
uqy(2,1)=0.01620327280964095
uqx(3,1)=-0.001914903012504294
uqy(3,1)=0.001122693458825686
uqx(3,2)=-0.001911747148401908
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0001755820633807187
uqy(4,1)=0.0001029603397676965
uqx(4,2)=-0.0002177119608478438
uqy(4,2)=5.372503557832448e-05
uqx(4,3)=-0.0002177119608478438
uqy(4,3)=-5.372503557832448e-05
uqx(5,1)=-1.860418596133933e-05
uqy(5,1)=1.090765125005818e-05
uqx(5,2)=-2.78198841681919e-05
uqy(5,2)=9.737957708444016e-06
uqx(5,3)=-2.811427798376512e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-2.78198841681919e-05
uqy(5,4)=-9.737957708444016e-06
endif
if(ntheta.EQ.18)then
uqx(2,1)=-0.0298059688900068
uqy(2,1)=0.01742692140290434
uqx(3,1)=-0.002261402797319784
uqy(3,1)=0.001324791852622696
uqx(3,2)=-0.00218979975575966
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0002263307798910572
uqy(4,1)=0.0001326134069799172
uqx(4,2)=-0.0002723167306780222
uqy(4,2)=6.749853746454043e-05
uqx(4,3)=-0.0002723167306780222
uqy(4,3)=-6.749853746454043e-05
uqx(5,1)=-2.618344775700557e-05
uqy(5,1)=1.533904799393413e-05
uqx(5,2)=-3.799806596259685e-05
uqy(5,2)=1.333194246897647e-05
uqx(5,3)=-3.804957109168239e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-3.799806596259685e-05
uqy(5,4)=-1.333194246897647e-05
endif
if(ntheta.EQ.19)then
uqx(2,1)=-0.03184839354005421
uqy(2,1)=0.01860686639767005
uqx(3,1)=-0.002632500871417905
uqy(3,1)=0.001541115305242353
uqx(3,2)=-0.00247086793286604
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0002854464623341651
uqy(4,1)=0.0001671332222719948
uqx(4,2)=-0.0003332456229776692
uqy(4,2)=8.294893631494263e-05
uqx(4,3)=-0.0003332456229776692
uqy(4,3)=-8.294893631494263e-05
uqx(5,1)=-3.5777331916109e-05
uqy(5,1)=2.094453998211611e-05
uqx(5,2)=-5.041053364039542e-05
uqy(5,2)=1.772747608189362e-05
uqx(5,3)=-5.003665347329104e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-5.041053364039542e-05
uqy(5,4)=-1.772747608189362e-05
endif
if(ntheta.EQ.20)then
uqx(2,1)=-0.03380921098354144
uqy(2,1)=0.01973915809034029
uqx(3,1)=-0.003025604178510975
uqy(3,1)=0.001770152796829164
uqx(3,2)=-0.002751053010400932
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0003530806004474008
uqy(4,1)=0.0002066047359070089
uqx(4,2)=-0.0003999707148415089
uqy(4,2)=9.995079644234232e-05
uqx(4,3)=-0.0003999707148415089
uqy(4,3)=-9.995079644234232e-05
uqx(5,1)=-4.761749968378625e-05
uqy(5,1)=2.785831220463549e-05
uqx(5,2)=-6.517066822235473e-05
uqy(5,2)=2.296892009302885e-05
uqx(5,3)=-6.414716777612842e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-6.517066822235473e-05
uqy(5,4)=-2.296892009302885e-05
endif
if(ntheta.EQ.21)then
uqx(2,1)=-0.03568421614472853
uqy(2,1)=0.0208215014107083
uqx(3,1)=-0.003438105278930305
uqy(3,1)=0.002010387168015211
uqx(3,2)=-0.003026890903307585
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0004292590682540321
uqy(4,1)=0.0002510400416038367
uqx(4,2)=-0.0004718571203498731
uqy(4,2)=0.0001183459294007554
uqx(4,3)=-0.0004718571203498731
uqy(4,3)=-0.0001183459294007554
uqx(5,1)=-6.190907493630538e-05
uqy(5,1)=3.61990173774744e-05
uqx(5,2)=-8.233825903866238e-05
uqy(5,2)=2.908155778293243e-05
uqx(5,3)=-8.040403158689358e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-8.233825903866238e-05
uqy(5,4)=-2.908155778293243e-05
endif
if(ntheta.EQ.22)then
uqx(2,1)=-0.03747156089826429
uqy(2,1)=0.02185294737277631
uqx(3,1)=-0.003867466589799954
uqy(3,1)=0.002260343560684938
uqx(3,2)=-0.003295398510176307
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0005139023086210218
uqy(4,1)=0.0003003901861872707
uqx(4,2)=-0.0005482018340347233
uqy(4,2)=0.0001379531343067409
uqx(4,3)=-0.0005482018340347233
uqy(4,3)=-0.0001379531343067409
uqx(5,1)=-7.882756813514927e-05
uqy(5,1)=4.606802067573352e-05
uqx(5,2)=-0.0001019216358938264
uqy(5,2)=3.607207962412092e-05
uqx(5,3)=-9.878582619642953e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0001019216358938264
uqy(5,4)=-3.607207962412092e-05
endif
if(ntheta.EQ.23)then
uqx(2,1)=-0.03917124800252712
uqy(2,1)=0.02283360020288869
uqx(3,1)=-0.00431127715331047
uqy(3,1)=0.002518621941051702
uqx(3,2)=-0.003554081631682294
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0006068462026936415
uqy(4,1)=0.0003545573425025866
uqx(4,2)=-0.0006282673934951395
uqy(4,2)=0.00015857707553736
uqx(4,3)=-0.0006282673934951395
uqy(4,3)=-0.00015857707553736
uqx(5,1)=-9.851794930567258e-05
uqy(5,1)=5.754889134915976e-05
uqx(5,2)=-0.0001238821974354894
uqy(5,2)=4.392995908857313e-05
uqx(5,3)=-0.0001192328797561733
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0001238821974354894
uqy(5,4)=-4.392995908857313e-05
endif
if(ntheta.EQ.24)then
uqx(2,1)=-0.04078467983132041
uqy(2,1)=0.02376435704408403
uqx(3,1)=-0.004767286855073334
uqy(3,1)=0.00278391658984503
uqx(3,2)=-0.003800916556729175
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0007078617681184562
uqy(4,1)=0.0004134062745389956
uqx(4,2)=-0.0007113093010572949
uqy(4,2)=0.0001800159036691447
uqx(4,3)=-0.0007113093010572949
uqy(4,3)=-0.0001800159036691447
uqx(5,1)=-0.0001210953749828021
uqy(5,1)=7.070785140253987e-05
uqx(5,2)=-0.0001481403717251601
uqy(5,2)=5.26293883724199e-05
uqx(5,3)=-0.0001416540841691958
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0001481403717251601
uqy(5,4)=-5.26293883724199e-05
endif
if(ntheta.EQ.25)then
uqx(2,1)=-0.04231427495534361
uqy(2,1)=0.02464668743107288
uqx(3,1)=-0.005233423406853966
uqy(3,1)=0.003055025655109403
uqx(3,2)=-0.004034315476632931
uqy(3,2)=0.0_dp
uqx(4,1)=-0.000816672621751051
uqy(4,1)=0.0004767744840141218
uqx(4,2)=-0.0007965971920023176
uqy(4,2)=0.0002020674940364024
uqx(4,3)=-0.0007965971920023176
uqy(4,3)=-0.0002020674940364024
uqx(5,1)=-0.0001466470646981058
uqy(5,1)=8.559488727257606e-05
uqx(5,2)=-0.0001745822216037299
uqy(5,2)=6.213149817457247e-05
uqx(5,3)=-0.0001659337332567695
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0001745822216037299
uqy(5,4)=-6.213149817457247e-05
endif
if(ntheta.EQ.26)then
uqx(2,1)=-0.0437631542688065
uqy(2,1)=0.02548245325111225
uqx(3,1)=-0.005707796948781806
uqy(3,1)=0.003330853585610454
uqx(3,2)=-0.00425308343192316
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0009329697453595618
uqy(4,1)=0.000544480777040355
uqx(4,2)=-0.0008834303440076407
uqy(4,2)=0.0002245343486965423
uqx(4,3)=-0.0008834303440076407
uqy(4,3)=-0.0002245343486965423
uqx(5,1)=-0.0001752348810775276
uqy(5,1)=0.0001022452652849156
uqx(5,2)=-0.0002030661159222708
uqy(5,2)=7.23866540137963e-05
uqx(5,3)=-0.0001919379139721673
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0002030661159222708
uqy(5,4)=-7.23866540137963e-05
endif
if(ntheta.EQ.27)then
uqx(2,1)=-0.04513489169976175
uqy(2,1)=0.02627376613853282
uqx(3,1)=-0.006188696292743292
uqy(3,1)=0.003610408772998517
uqx(3,2)=-0.0044563721863212
uqy(3,2)=0.0_dp
uqx(4,1)=-0.001056423499571107
uqy(4,1)=0.0006163322244658482
uqx(4,2)=-0.0009711483975435052
uqy(4,2)=0.0002472273003311512
uqx(4,3)=-0.0009711483975435052
uqy(4,3)=-0.0002472273003311512
uqx(5,1)=-0.0002068982555902986
uqy(5,1)=0.0001206812433280477
uqx(5,2)=-0.0002334290810502318
uqy(5,2)=8.333668679649319e-05
uqx(5,3)=-0.0002195201853748941
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0002334290810502318
uqy(5,4)=-8.333668679649319e-05
endif
if(ntheta.EQ.28)then
uqx(2,1)=-0.04643332085491404
uqy(2,1)=0.02702287718328207
uqx(3,1)=-0.006674579953450168
uqy(3,1)=0.003892798220305804
uqx(3,2)=-0.004643634524000381
uqy(3,2)=0.0_dp
uqx(4,1)=-0.001186693080114393
uqy(4,1)=0.0006921296326732865
uqx(4,2)=-0.001059138218690306
uqy(4,2)=0.0002699681973246943
uqx(4,3)=-0.001059138218690306
uqy(4,3)=-0.0002699681973246943
uqx(5,1)=-0.0002416571953113832
uqy(5,1)=0.000140913825684749
uqx(5,2)=-0.0002654926052230902
uqy(5,2)=9.491697014352015e-05
uqx(5,3)=-0.000248526431033852
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0002654926052230902
uqy(5,4)=-9.491697014352015e-05
endif
if(ntheta.EQ.29)then
uqx(2,1)=-0.0476623722906518
uqy(2,1)=0.02773208419857225
uqx(3,1)=-0.007164065345384643
uqy(3,1)=0.004177221193478502
uqx(3,2)=-0.004814580805092504
uqy(3,2)=0.0_dp
uqx(4,1)=-0.001323433931627007
uqy(4,1)=0.0007716718280408433
uqx(4,2)=-0.001146837865060967
uqy(4,2)=0.0002925918557552575
uqx(4,3)=-0.001146837865060967
uqy(4,3)=-0.0002925918557552575
uqx(5,1)=-0.0002795152299791238
uqy(5,1)=0.0001629444802374371
uqx(5,2)=-0.0002990678227357772
uqy(5,2)=0.0001070583234817996
uqx(5,3)=-0.0002787988860265953
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0002990678227357772
uqy(5,4)=-0.0001070583234817996
endif
if(ntheta.EQ.30)then
uqx(2,1)=-0.04882563370160641
uqy(2,1)=0.02840347765698087
uqx(3,1)=-0.00765593696089868
uqy(3,1)=0.004462974022425572
uqx(3,2)=-0.004969132784832358
uqy(3,2)=0.0_dp
uqx(4,1)=-0.00146630725510065
uqy(4,1)=0.0008547611790123954
uqx(4,2)=-0.001233739988517768
uqy(4,2)=0.0003149494088037902
uqx(4,3)=-0.001233739988517768
uqy(4,3)=-0.0003149494088037902
uqx(5,1)=-0.0003204630115092076
uqy(5,1)=0.0001867672361436166
uqx(5,2)=-0.0003339607774787402
uqy(5,2)=0.0001196892313188401
uqx(5,3)=-0.0003101796338676208
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0003339607774787402
uqy(5,4)=-0.0001196892313188401
endif
if(ntheta.EQ.31)then
uqx(2,1)=-0.04992193062183704
uqy(2,1)=0.0290363736961492
uqx(3,1)=-0.008149400568775942
uqy(3,1)=0.004749598900816575
uqx(3,2)=-0.005107282626575934
uqy(3,2)=0.0_dp
uqx(4,1)=-0.001615038562818766
uqy(4,1)=0.0009412378657909575
uqx(4,2)=-0.001319412580829482
uqy(4,2)=0.0003369375238227954
uqx(4,3)=-0.001319412580829482
uqy(4,3)=-0.0003369375238227954
uqx(5,1)=-0.0003644932065908745
uqy(5,1)=0.0002123773891528697
uqx(5,2)=-0.0003699867939099062
uqy(5,2)=0.0001327448853399043
uqx(5,3)=-0.0003425164020803855
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0003699867939099062
uqy(5,4)=-0.0001327448853399043
endif
if(ntheta.EQ.32)then
uqx(2,1)=-0.05090329799667861
uqy(2,1)=0.02960289728720748
uqx(3,1)=-0.008646542228826316
uqy(3,1)=0.00503832273936899
uqx(3,2)=-0.005227869327712846
uqy(3,2)=0.0_dp
uqx(4,1)=-0.001769947347072314
uqy(4,1)=0.001031290132427384
uqx(4,2)=-0.001403673062984919
uqy(4,2)=0.0003587752109830483
uqx(4,3)=-0.001403673062984919
uqy(4,3)=-0.0003587752109830483
uqx(5,1)=-0.0004117279745827106
uqy(5,1)=0.0002398460686580247
uqx(5,2)=-0.0004070748751161033
uqy(5,2)=0.0001462439658497053
uqx(5,3)=-0.0003756924781756205
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0004070748751161033
uqy(5,4)=-0.0001462439658497053
endif
if(ntheta.EQ.33)then
uqx(2,1)=-0.0513606096830928
uqy(2,1)=0.02986539511944403
uqx(3,1)=-0.009170005042813718
uqy(3,1)=0.005342398460657084
uqx(3,2)=-0.005318370702859235
uqy(3,2)=0.0_dp
uqx(4,1)=-0.001935879097096427
uqy(4,1)=0.001127757228004547
uqx(4,2)=-0.001487689894477332
uqy(4,2)=0.0003830788831587224
uqx(4,3)=-0.001487689894477332
uqy(4,3)=-0.0003830788831587224
uqx(5,1)=-0.000463394046633049
uqy(5,1)=0.0002698905192682427
uqx(5,2)=-0.0004460132086709931
uqy(5,2)=0.0001608687247518319
uqx(5,3)=-0.0004098000837104587
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0004460132086709931
uqy(5,4)=-0.0001608687247518319
endif
if(ntheta.EQ.34)then
uqx(2,1)=-0.04874144182744936
uqy(2,1)=0.0283417826549473
uqx(3,1)=-0.009858487467627138
uqy(3,1)=0.005742969922777688
uqx(3,2)=-0.005287936990165755
uqy(3,2)=0.0_dp
uqx(4,1)=-0.002144796524776807
uqy(4,1)=0.001249360187141939
uqx(4,2)=-0.001578285629772355
uqy(4,2)=0.0004267732968133611
uqx(4,3)=-0.001578285629772355
uqy(4,3)=-0.0004267732968133611
uqx(5,1)=-0.0005276562550504655
uqy(5,1)=0.000307285033763345
uqx(5,2)=-0.0004925301179564875
uqy(5,2)=0.0001813270867116877
uqx(5,3)=-0.0004456179148133969
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0004925301179564875
uqy(5,4)=-0.0001813270867116877
endif
if(ntheta.EQ.35)then
uqx(2,1)=-0.03297296564955147
uqy(2,1)=0.01918709821552975
uqx(3,1)=-0.01123505415890843
uqy(3,1)=0.006546108167901785
uqx(3,2)=-0.004684706983380724
uqy(3,2)=0.0_dp
uqx(4,1)=-0.002533504540963353
uqy(4,1)=0.001476247589033027
uqx(4,2)=-0.001676066560316715
uqy(4,2)=0.0005604350432240707
uqx(4,3)=-0.001676066560316715
uqy(4,3)=-0.0005604350432240707
uqx(5,1)=-0.0006414180926895605
uqy(5,1)=0.0003736003362540499
uqx(5,2)=-0.0005690638663033982
uqy(5,2)=0.0002272281056786763
uqx(5,3)=-0.0004791697292710052
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0005690638663033982
uqy(5,4)=-0.0002272281056786763
endif

!!!!!!!!!!!!!!!!!!! Layer 1 at z=0

 do i=1,ndim/2

   rx=coord(i,1)
   ry=coord(i,2)

   delta(1)=(0.d0,0.d0)
   delta(2)=(0.d0,0.d0)

!write(*,*) n,m,vk(1),vk(2),rx,ry

   ncount=0
     do nn=1,5,2
       do n=2,5
         do m=1,n-1

           vk(1)=(n-1)*Gn(nn,1)+(m-1)*Gn(nn+2,1)
           vk(2)=(n-1)*Gn(nn,2)+(m-1)*Gn(nn+2,2)

           ux=cos((nn-1)*pi/3.0_dp)*uqx(n,m)-sin((nn-1)*pi/3.0_dp)*uqy(n,m)
           uy=sin((nn-1)*pi/3.0_dp)*uqx(n,m)+cos((nn-1)*pi/3.0_dp)*uqy(n,m)

           delta(1)=delta(1)-2*ux*sin(vk(1)*rx+vk(2)*ry)
           delta(2)=delta(2)-2*uy*sin(vk(1)*rx+vk(2)*ry)

         end do
       end do
     end do

!     if(coord(i,3).LT.deltaR) then
       coord(i,1)=rx+delta(1)/2
       coord(i,2)=ry+delta(2)/2
!     else


    coord(i+ndim/2,1)=-coord(i,1)!rx-delta(1)/2
    coord(i+ndim/2,2)=coord(i,2) !ry-delta(2)/2


enddo


end subroutine Relaxation_2



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine samplePoints_dB(npointsBZ,npointsBack,numkPlot,ncount)

  integer(dp), intent(in) :: numkPlot
  integer(dp), intent(in) :: ncount
  integer(dp), intent(inout) :: npointsBZ(ncount,2),npointsBack(numkPlot+1,numkPlot+1)

  integer(dp) :: ivk,ivk1,ivk2,icount

  icount=0
  do ivk1=0,numkPlot
  do ivk2=0,numkPlot
    icount=icount+1
    npointsBZ(icount,1)=ivk1
    npointsBZ(icount,2)=ivk2
    npointsBack(npointsBZ(icount,1)+1,npointsBZ(icount,2)+1)=icount
  enddo
  enddo


  return
end subroutine samplePoints_dB


subroutine samplePoints1(npointsBZ,npointsBack,numk,ncount)

  integer(dp), intent(in) :: numk
  integer(dp), intent(in) :: ncount
  integer(dp), intent(inout) :: npointsBZ(ncount,8),npointsBack(numk,numk)

  integer(dp) :: ivk,ivk1,ivk2,icount

  icount=0
  do ivk1=0,numk-1
  do ivk2=0,numk-1
    icount=icount+1
    npointsBZ(icount,1)=ivk1
    npointsBZ(icount,2)=ivk2
    npointsBZ(icount,3)=0
    npointsBZ(icount,4)=1
    npointsBack(npointsBZ(icount,1)+1,npointsBZ(icount,2)+1)=icount
  enddo
  enddo


  return
end subroutine samplePoints1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module ResponseGeometry
