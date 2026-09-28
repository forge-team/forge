! Adapted from inbox/Hamiltonian.f90. Physical settings come from ../Setup.f90.
! Screening is shared with the solver; observable formulas are retained.
 module Hamiltonian

 use lapack_routines
 use HartreeFock, only: fv, LongRangeInteraction
 use Setup, only: dp, numk, ntheta, ndim, numI, numC, numS, nfilling, nrelax, xi, fphase, pressure, U, a0, pi, tz, a1, a2, dirFock

 implicit none
 !Fock-SYM-filling+0.0-i22-relax1-eps4.0-U3.10-scr1-xi100.0-delta.000-numI2-numk6-dp8.dat
 integer(dp), parameter :: numb = ndim
 integer(dp), parameter :: ncb = ndim/2
 integer(dp), parameter :: nFermiSea=1!ndim/2-numb+ncb+1-128!1<=nFermiSea=<ndim/2-numb+ncb+1
 integer(dp), parameter :: nRandom=1
 integer(dp), parameter :: nDrude=0_dp
 integer(dp), parameter :: nDos=0_dp
 !integer(dp) , parameter :: nT=7
 !real(dp)    , parameter :: ftemperature(nT) =[0._dp,1._dp,10._dp,36._dp,50._dp,200._dp,300._dp]
 integer(dp) , parameter :: nT=1
 real(dp)    , parameter :: ftemperature(nT) =[1._dp]
 integer(dp), parameter :: nmu=40

 integer(dp), parameter :: nPrintDensity=0
 integer(dp), parameter :: nRead=1
 integer(dp), parameter :: nWrite=1
 character(8) :: dir='output4/'
 real(dp), parameter :: deltaR = 0.1_dp
 real(dp), parameter :: r0 = 0.184_dp
 real(dp), parameter :: d0 = 1.35772_dp/pressure
 real(dp), parameter :: tperp = 1._dp
 real(dp), parameter :: mass1=0._dp!-0.005_dp

 real(dp), parameter :: T_in_ev = 8.617342791d-5
 real(dp), parameter :: sq3 = sqrt(3.0_dp)
 real(dp), parameter :: vkd=4.0_dp*pi/3.0_dp


 !real(dp) :: vGamma(2),vK1(2),vK2(2),vM(2),vshift(2)
 !real(dp) :: aGM,aKG,aMK,aT,eD

contains


pure elemental function fermi_dist(e,t)

 real(dp), intent(in) :: e
 real(dp), optional, intent(in) :: t
 real(dp) :: fermi_dist

 if(present(t) .and. t > 1d-10) then
   fermi_dist = 1.0_dp/( 1.0_dp + exp(e/(t*t_in_ev)) )
 else
   if( e<0.0_dp ) then
     fermi_dist = 1.0_dp
   else
     fermi_dist = 0.0_dp
   end if
 end if

end function fermi_dist

!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!


subroutine longRange(fLongRange,coord,ndim,t1,t2)
integer(dp), intent(in) :: ndim
real(dp), intent(in) :: coord(ndim,3), t1(2), t2(2)
real(dp), intent(inout) :: fLongRange(ndim,ndim)
call LongRangeInteraction(fLongRange,coord,ndim,t1,t2)
end subroutine longRange
!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!


!c!c!cC
!c!c!cC
!c!c!cC

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

function sigma_func(z1,z2,energy)

 real(dp)   , intent(in) :: energy
 complex(dp), intent(in) :: z1, z2
 real(dp) :: sigma_func

 sigma_func = real(conjg(z1)*z2,dp)*energy

end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function temperature(n)

 integer(dp), intent(in) :: n
 real(dp) :: temperature

 if(n==1_dp) then
   temperature = 0.0000861733_dp
 else
   temperature = 300.0_dp*0.0000861733_dp
 end if

end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine ComplexH0(zh,coord,ndim,vk,tn)

 integer(dp), intent(in)    :: ndim
 real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3)
 complex(dp), intent(out) :: zh(ndim,ndim)

 integer(dp) :: i, j, nt
 real(dp)    :: rij, coord_diff
 complex(dp) :: zi_vk_tn,zphase


zh = cmplx(0.0_dp,0.0_dp,dp)

 do i = 1 , ndim
   do j = 1 , i

     zphase=exp(cmplx(0.0_dp,-dot_product(vk,coord(i,1:2)-coord(j,1:2))*fphase,dp))

     rij = norm2(coord(i,:)-coord(j,:))

     coord_diff = (coord(i,3)-coord(j,3))**2

     if(coord_diff < deltaR .and. rij > deltaR) then
       zh(i,j) = ft(rij)*zphase
     endif

     if(coord_diff > deltaR) then
        zh(i,j) = ftperp(rij)*zphase
     endif

!!!!!! Coupling to adjacent Wigner-Seitz cell with tn
     do nt = 1 , 3
       rij = norm2(coord(i,:)-[coord(j,1:2)-tn(:,nt),coord(j,3)])
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
       if(coord_diff < deltaR)then
         zh(i,j)=zh(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase
       else
         zh(i,j)=zh(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase
       endif
!!!!!! Coupling to adjacent Wigner-Seitz cell with -tn
  rij = norm2(coord(i,:)-[coord(j,1:2)+tn(:,nt),coord(j,3)])
  zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,tn(:,nt)),dp)
  if(coord_diff < deltaR)then
    zh(i,j)=zh(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase
  else
    zh(i,j)=zh(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase
  endif
enddo

   end do
 end do

 do i=1,ndim/4
 zh(i,i)=zh(i,i)+mass1
 j=i+ndim/4
 zh(j,j)=zh(j,j)-mass1
 enddo

end subroutine ComplexH0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine ComplexH0_phase(zh,coord,ind_j,ind_l,ndim,ind,vk,tn)

 integer(dp), intent(in)    :: ndim,ind
 integer(dp), intent(in)    :: ind_j(ind),ind_l(ind)
 real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3)
 complex(dp), intent(out) :: zh(ndim,ndim)

 integer(dp) :: i, j, n1,n2,icount,ncount
 real(dp)    :: rij, coord_diff,coord2(ndim,2)
 complex(dp) :: zi_vk_tn,zphase


zh = cmplx(0.0_dp,0.0_dp,dp)

coord2(:,1:2)=coord(:,1:2)

ncount=0

 do i = 1 , ndim
   do j = 1 , i

     zphase=exp(cmplx(0.0_dp,-dot_product(vk,coord2(i,:)-coord2(j,:))*fphase,dp))

     rij = norm2(coord(i,:)-coord(j,:))

     coord_diff = (coord(i,3)-coord(j,3))**2

     if(coord_diff < deltaR .and. rij > deltaR) then
       zh(i,j) = ft(rij)*zphase
     endif

     if(coord_diff > deltaR) then
        zh(i,j) = ftperp(rij)*zphase
     endif

     do icount=2,ind
     n1=ind_j(icount)
     n2=ind_l(icount)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
   if(coord_diff < deltaR)then
      zh(i,j)=zh(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase
   else
      zh(i,j)=zh(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase
   endif

   rij = norm2(coord(i,:)-[coord(j,1:2)+n1*tn(:,1)+n2*tn(:,2),coord(j,3)])

   zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
   if(coord_diff < deltaR)then
      zh(i,j)=zh(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase
   else
      zh(i,j)=zh(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase
   endif

   enddo

   end do
 end do

 !write(*,*) ncount/real(ndim*ndim)

 do i=1,ndim/4
 zh(i,i)=zh(i,i)+mass1
 j=i+ndim/4
 zh(j,j)=zh(j,j)-mass1
 enddo


end subroutine ComplexH0_phase

!c!c!cC
!c!c!cC
!c!c!cC

subroutine ComplexH0_MF(zh,coord,potential,alpha,zfock,ndim)

integer(dp), intent(in)    :: ndim

real(dp)   , intent(in)    :: coord(ndim,3), potential(ndim),alpha
complex(dp) , intent(in)   :: zfock(ndim,ndim)
complex(dp), intent(out) :: zh(ndim,ndim)

integer(dp) :: i, j, nt
real(dp)    :: rij, coord_diff
complex(dp) :: zi_vk_tn,zphase


zh(:,:) = cmplx(0.0_dp,0.0_dp,dp)

 do i = 1 , ndim
   zh(i,i)=potential(i)
   do j = 1 , i

   rij = norm2(coord(i,:)-coord(j,:))

   coord_diff = (coord(i,3)-coord(j,3))**2

   zh(i,j)=zh(i,j)-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j))


   if(coord_diff < deltaR .and. rij > deltaR) then
     zh(i,j) = zh(i,j)+ft(rij)
   endif

   if(coord_diff > deltaR) then
      zh(i,j) = zh(i,j) + ftperp(rij)
   endif

   end do
 end do


end subroutine ComplexH0_MF

!c!c!cC
!c!c!cC
!c!c!cC

subroutine ComplexHk_MF_UnitCell(zh,coord,potential,alpha,zfock,ind_j,ind_l,ndim,ind,vk,tn)

integer(dp), intent(in)    :: ndim, ind

integer(dp), intent(in)    :: ind_j(ind), ind_l(ind)
real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3),potential(ndim),alpha
complex(dp) , intent(in)   :: zfock(ndim,ndim,ind)
complex(dp), intent(out) :: zh(ndim,ndim)

integer(dp) :: i, j, nt, icount, n1, n2
real(dp)    :: rij, coord_diff
complex(dp) :: zi_vk_tn,zsum,zphase

do i=1,ind
! write(*,*) i,vk(1),vk(2)
! write(*,*) i,zfock(345,140,i)
! write(*,*) i,conjg(zfock(345,140,i))
enddo


zh(:,:) = cmplx(0.0_dp,0.0_dp,dp)

 do i = 1 , ndim
   zh(i,i)=zh(i,i)+potential(i)
   do j = 1 , i

     zphase=exp(cmplx(0.0_dp,-dot_product(vk,coord(i,1:2)-coord(j,1:2))*fphase,dp))

rij = norm2(coord(i,:)-coord(j,:))

coord_diff = (coord(i,3)-coord(j,3))**2

zh(i,j)=zh(i,j)-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase


if(coord_diff < deltaR .and. rij > deltaR) then
  zh(i,j) = zh(i,j)+ft(rij)*zphase
endif

if(coord_diff > deltaR) then
   zh(i,j) = zh(i,j) + ftperp(rij)*zphase
endif


!!!!!! Coupling to adjacent Wigner-Seitz cell with tn
     do nt = 1 , 3
       rij = norm2(coord(i,:)-[coord(j,1:2)-tn(:,nt),coord(j,3)])

       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)))

       if(coord_diff < deltaR)then
         zh(i,j)=zh(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase
       else
         zh(i,j)=zh(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase
       endif

!!!!!! Coupling to adjacent Wigner-Seitz cell with tn

  rij = norm2(coord(i,:)-[coord(j,1:2)+tn(:,nt),coord(j,3)])

  zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,tn(:,nt)),dp)

  if(coord_diff < deltaR)then
    zh(i,j)=zh(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase
  else
    zh(i,j)=zh(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase
  endif

     end do


   zsum=cmplx(0.0_dp,0.0_dp,dp)
   do icount=2,ind
   n1=ind_j(icount)
   n2=ind_l(icount)

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])

   zsum=zsum-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)

   rij = norm2(coord(j,:)-[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])

   zsum=zsum-alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)

   enddo

   zh(i,j)=zh(i,j)+zsum*zphase

   end do
 end do

 do i=1,ndim/4
 zh(i,i)=zh(i,i)+mass1
 j=i+ndim/4
 zh(j,j)=zh(j,j)-mass1
 enddo


end subroutine ComplexHk_MF_UnitCell

!c!c!cC
!c!c!cC
!c!c!cC

subroutine ComplexHk_MF(zh,coord,potential,alpha,zfock,ind_j,ind_l,ndim,ind,vk,tn)

integer(dp), intent(in)    :: ndim, ind

integer(dp), intent(in)    :: ind_j(ind), ind_l(ind)
real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3),potential(ndim),alpha
complex(dp) , intent(in)   :: zfock(ndim,ndim,ind)
complex(dp), intent(out) :: zh(ndim,ndim)

integer(dp) :: i, j, nt, icount, n1, n2
real(dp)    :: rij, coord_diff
complex(dp) :: zi_vk_tn,zsum,zphase

do i=1,ind
! write(*,*) i,vk(1),vk(2)
! write(*,*) i,zfock(345,140,i)
! write(*,*) i,conjg(zfock(345,140,i))
enddo


zh(:,:) = cmplx(0.0_dp,0.0_dp,dp)

 do i = 1 , ndim
   zh(i,i)=potential(i)
   do j = 1 , i

     zphase=exp(cmplx(0.0_dp,-dot_product(vk,coord(i,1:2)-coord(j,1:2))*fphase,dp))

rij = norm2(coord(i,:)-coord(j,:))

coord_diff = (coord(i,3)-coord(j,3))**2

zh(i,j)=zh(i,j)-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase


if(coord_diff < deltaR .and. rij > deltaR) then
  zh(i,j) = zh(i,j)+ft(rij)*zphase
endif

if(coord_diff > deltaR) then
   zh(i,j) = zh(i,j) + ftperp(rij)*zphase
endif

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

     do icount=2,ind
     n1=ind_j(icount)
     n2=ind_l(icount)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
   if(coord_diff < deltaR)then
      zh(i,j)=zh(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase
   else
      zh(i,j)=zh(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase
   endif

   rij = norm2(coord(i,:)-[coord(j,1:2)+n1*tn(:,1)+n2*tn(:,2),coord(j,3)])

   zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
   if(coord_diff < deltaR)then
      zh(i,j)=zh(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase
   else
      zh(i,j)=zh(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase
   endif

   enddo

!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!

   zsum=cmplx(0.0_dp,0.0_dp,dp)
   do icount=2,ind
   n1=ind_j(icount)
   n2=ind_l(icount)

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])

   zsum=zsum-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)

   rij = norm2(coord(j,:)-[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])

   zsum=zsum-alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)

   enddo

   zh(i,j)=zh(i,j)+zsum*zphase

   end do
 end do

 !do i=1,ndim/4
 !zh(i,i)=zh(i,i)+mass1
 !j=i+ndim/4
 !zh(j,j)=zh(j,j)-mass1
 !ssenddo



end subroutine ComplexHk_MF

!c!c!cC
!c!c!cC
!c!c!cC

subroutine vxy_MF(zvxy,coord,alpha,zfock,ind_j,ind_l,ndim,ind,vk,tn,nxy)

integer(dp), intent(in)    :: ndim, ind, nxy

integer(dp), intent(in)    :: ind_j(ind), ind_l(ind)
real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3),alpha
complex(dp) , intent(in)   :: zfock(ndim,ndim,ind)
complex(dp), intent(out) :: zvxy(ndim,ndim)

integer(dp) :: i, j, nt, icount, n1, n2
real(dp)    :: rij, rxy, coord_diff
complex(dp) :: zi_vk_tn,zsum,zphase


 zvxy = cmplx(0.0_dp,0.0_dp,dp)

 do i = 1 , ndim
   do j = i , ndim

     zphase=exp(cmplx(0.0_dp,-dot_product(vk,coord(i,1:2)-coord(j,1:2))*fphase,dp))

     rij = norm2(coord(i,:)-coord(j,:))
     rxy = coord(i,nxy)-coord(j,nxy)

     coord_diff = (coord(i,3)-coord(j,3))**2

     zvxy(i,j)=zvxy(i,j)-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*cmplx(0.0_dp,rxy,dp)


if(coord_diff < deltaR .and. rij > deltaR) then
  zvxy(i,j)=zvxy(i,j)+ft(rij)*zphase*cmplx(0.0_dp,rxy,dp)
endif

if(coord_diff > deltaR) then
   zvxy(i,j)=zvxy(i,j) + ftperp(rij)*zphase*cmplx(0.0_dp,rxy,dp)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

     do icount=2,ind
     n1=ind_j(icount)
     n2=ind_l(icount)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)
   !rxy = coord(i,nxy)-coord(j,nxy)

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
   if(coord_diff < deltaR)then
      zvxy(i,j)=zvxy(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   else
      zvxy(i,j)=zvxy(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   endif

   rij = norm2(coord(i,:)-[coord(j,1:2)+n1*tn(:,1)+n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)-n1*tn(nxy,1)-n2*tn(nxy,2)

   zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
   if(coord_diff < deltaR)then
      zvxy(i,j)=zvxy(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   else
      zvxy(i,j)=zvxy(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   endif

   enddo

!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!

   zsum=cmplx(0.0_dp,0.0_dp,dp)
   do icount=2,ind
   n1=ind_j(icount)
   n2=ind_l(icount)

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   zsum=zsum-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)

   rij = norm2(coord(j,:)-[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])
   rxy=coord(j,nxy)-coord(i,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   zsum=zsum-alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*cmplx(0.0_dp,rxy,dp)

   enddo

   zvxy(i,j)=zvxy(i,j)+zsum*zphase

   end do
 end do



end subroutine vxy_MF

!c!c!cC
!c!c!cC
!c!c!cC

subroutine vxy12I_MF(zvxy1,zvxy2,zvxyI,coord,alpha,zfock,ind_j,ind_l,ndim,ind,vk,tn,nxy)

integer(dp), intent(in)    :: ndim, ind, nxy

integer(dp), intent(in)    :: ind_j(ind), ind_l(ind)
real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3),alpha
complex(dp) , intent(in)   :: zfock(ndim,ndim,ind)
complex(dp), intent(out) :: zvxy1(:,:),zvxy2(:,:),zvxyI(:,:)

integer(dp) :: i, j, nt, icount, n1, n2
real(dp)    :: rij, rxy, coord_diff
complex(dp) :: zi_vk_tn,zsum1,zsum2,zsumI,zphase


 zvxy1 = cmplx(0.0_dp,0.0_dp,dp)
 zvxy2 = cmplx(0.0_dp,0.0_dp,dp)
 zvxyI = cmplx(0.0_dp,0.0_dp,dp)


 do i = 1 , ndim
   do j = i , ndim

     zphase=exp(cmplx(0.0_dp,-dot_product(vk,coord(i,1:2)-coord(j,1:2))*fphase,dp))

     rij = norm2(coord(i,:)-coord(j,:))
     rxy = coord(i,nxy)-coord(j,nxy)


if(coord(i,3)<deltaR .and. coord(j,3)<deltaR .and. rij>deltaR) then
  zvxy1(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*cmplx(0.0_dp,rxy,dp)
  zvxy1(i,j)=zvxy1(i,j)+ft(rij)*zphase*cmplx(0.0_dp,rxy,dp)
else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR .and. rij>deltaR) then
  zvxy2(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*cmplx(0.0_dp,rxy,dp)
  zvxy2(i,j)=zvxy2(i,j)+ft(rij)*zphase*cmplx(0.0_dp,rxy,dp)
else if(coord(i,3).NE.coord(j,3)) then
  zvxyI(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*cmplx(0.0_dp,rxy,dp)
  zvxyI(i,j)=zvxyI(i,j) + ftperp(rij)*zphase*cmplx(0.0_dp,rxy,dp)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

     do icount=2,ind
     n1=ind_j(icount)
     n2=ind_l(icount)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
      zvxy1(i,j)=zvxy1(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
      zvxy2(i,j)=zvxy2(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3).NE.coord(j,3)) then
      zvxyI(i,j)=zvxyI(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   endif

   rij = norm2(coord(i,:)-[coord(j,1:2)+n1*tn(:,1)+n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)-n1*tn(nxy,1)-n2*tn(nxy,2)

   zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
      zvxy1(i,j)=zvxy1(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
      zvxy2(i,j)=zvxy2(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3).NE.coord(j,3)) then
      zvxyI(i,j)=zvxyI(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   endif

   enddo

!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!

   zsum1=cmplx(0.0_dp,0.0_dp,dp)
   zsum2=cmplx(0.0_dp,0.0_dp,dp)
   zsumI=cmplx(0.0_dp,0.0_dp,dp)
   do icount=2,ind
   n1=ind_j(icount)
   n2=ind_l(icount)

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
      zsum1=zsum1-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
      zsum2=zsum2-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3).NE.coord(j,3)) then
      zsumI=zsumI-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
   endif

   rij = norm2(coord(j,:)-[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])
   rxy=coord(j,nxy)-coord(i,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
   zsum1=zsum1+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
   zsum2=zsum2+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3).NE.coord(j,3)) then
   zsumI=zsumI+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
   endif

   enddo

   zvxy1(i,j)=zvxy1(i,j)+zsum1*zphase
   zvxy2(i,j)=zvxy2(i,j)+zsum2*zphase
   zvxyI(i,j)=zvxyI(i,j)+zsumI*zphase

   end do
 end do



end subroutine vxy12I_MF

!c!c!cC
!c!c!cC
!c!c!cC

subroutine vxy12_MF(zvxy,coord,alpha,zfock,ind_j,ind_l,ndim,ind,vk,tn,nxy)

integer(dp), intent(in)    :: ndim, ind, nxy

integer(dp), intent(in)    :: ind_j(ind), ind_l(ind)
real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3),alpha
complex(dp) , intent(in)   :: zfock(ndim,ndim,ind)
complex(dp), intent(out) :: zvxy(:,:)

integer(dp) :: i, j, nt, icount, n1, n2
real(dp)    :: rij, rxy, coord_diff
complex(dp) :: zi_vk_tn,zsum1,zsum2,zsumI,zphase


 zvxy = cmplx(0.0_dp,0.0_dp,dp)


 do i = 1 , ndim
   do j = i , ndim

     zphase=exp(cmplx(0.0_dp,-dot_product(vk,coord(i,1:2)-coord(j,1:2))*fphase,dp))

     rij = norm2(coord(i,:)-coord(j,:))
     rxy = coord(i,nxy)-coord(j,nxy)


if(coord(i,3)<deltaR .and. coord(j,3)<deltaR .and. rij>deltaR) then
  zvxy(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*cmplx(0.0_dp,rxy,dp)
  zvxy(i,j)=zvxy(i,j)+ft(rij)*zphase*cmplx(0.0_dp,rxy,dp)
else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR .and. rij>deltaR) then
  zvxy(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*cmplx(0.0_dp,rxy,dp)
  zvxy(i,j)=zvxy(i,j)+ft(rij)*zphase*cmplx(0.0_dp,rxy,dp)
else if(coord(i,3).NE.coord(j,3)) then
  zvxy(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*cmplx(0.0_dp,rxy,dp)
  zvxy(i,j)=zvxy(i,j) + ftperp(rij)*zphase*cmplx(0.0_dp,rxy,dp)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

     do icount=2,ind
     n1=ind_j(icount)
     n2=ind_l(icount)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
      zvxy(i,j)=zvxy(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
      zvxy(i,j)=zvxy(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3).NE.coord(j,3)) then
      zvxy(i,j)=zvxy(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   endif

   rij = norm2(coord(i,:)-[coord(j,1:2)+n1*tn(:,1)+n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)-n1*tn(nxy,1)-n2*tn(nxy,2)

   zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
      zvxy(i,j)=zvxy(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
      zvxy(i,j)=zvxy(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3).NE.coord(j,3)) then
      zvxy(i,j)=zvxy(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase*cmplx(0.0_dp,rxy,dp)
   endif

   enddo

!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!

   zsum1=cmplx(0.0_dp,0.0_dp,dp)
   zsum2=cmplx(0.0_dp,0.0_dp,dp)
   zsumI=cmplx(0.0_dp,0.0_dp,dp)
   do icount=2,ind
   n1=ind_j(icount)
   n2=ind_l(icount)

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
      zsum1=zsum1-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
      zsum2=zsum2-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3).NE.coord(j,3)) then
      zsumI=zsumI-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
   endif

   rij = norm2(coord(j,:)-[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])
   rxy=coord(j,nxy)-coord(i,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
   zsum1=zsum1+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
   zsum2=zsum2+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
   else if(coord(i,3).NE.coord(j,3)) then
   zsumI=zsumI+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
   endif

   enddo

   zvxy(i,j)=zvxy(i,j)+zsum1*zphase+zsum2*zphase+zsumI*zphase

   end do
 end do



end subroutine vxy12_MF



!c!c!cC
!c!c!cC
!c!c!cC

subroutine ComplexH0_Hartree(zh,coord,potential,ndim)

integer(dp), intent(in)    :: ndim

real(dp)   , intent(in)    :: coord(ndim,3), potential(ndim)
complex(dp), intent(out) :: zh(ndim,ndim)

integer(dp) :: i, j, nt
real(dp)    :: rij, coord_diff
complex(dp) :: zi_vk_tn


zh(:,:) = cmplx(0.0_dp,0.0_dp,dp)

 do i = 1 , ndim
   zh(i,i)=potential(i)
   do j = 1 , ndim

   rij = norm2(coord(i,:)-coord(j,:))

   coord_diff = (coord(i,3)-coord(j,3))**2

   if(coord_diff < deltaR .and. rij > deltaR) then
     zh(i,j) = zh(i,j)+ft(rij)
   endif

   if(coord_diff > deltaR) then
      zh(i,j) = zh(i,j) + ftperp(rij)
   endif


   end do
 end do

end subroutine ComplexH0_Hartree

!c!c!cC
!c!c!cC
!c!c!cC

subroutine ComplexHk_Hartree(zh,coord,ndim,vk,tn)

integer(dp), intent(in)    :: ndim

real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3)
complex(dp), intent(inout) :: zh(ndim,ndim)

integer(dp) :: i, j, nt, icount, n1, n2
real(dp)    :: rij, coord_diff
complex(dp) :: zi_vk_tn


zh(:,:) = cmplx(0.0_dp,0.0_dp,dp)

 do i = 1 , ndim
   do j = 1 , ndim

     coord_diff = (coord(i,3)-coord(j,3))**2

!!!!!! Coupling to adjacent Wigner-Seitz cell with tn
     do nt = 1 , 3
       rij = norm2(coord(i,:)-[coord(j,1:2)-tn(:,nt),coord(j,3)])

       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)

       if(coord_diff < deltaR)then
         zh(i,j)=zh(i,j)+ft(rij)*exp( zi_vk_tn)
         zh(j,i)=zh(j,i)+ft(rij)*exp(-zi_vk_tn)
       else
         zh(i,j)=zh(i,j)+ftperp(rij)*exp( zi_vk_tn)
         zh(j,i)=zh(j,i)+ftperp(rij)*exp(-zi_vk_tn)
       endif

     end do

   end do
 end do

end subroutine ComplexHk_Hartree

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine ComplexH0_sym(zh,coord,ndim,vk,tn)

 integer(dp), intent(in)    :: ndim
 real(dp)   , intent(in)    :: tn(2,3), vk(2)
 real(dp)   , intent(in)    :: coord(ndim,3)
 complex(dp), intent(inout) :: zh(ndim,ndim)

 integer(dp) :: i, j, nt, l1, l2
 real(dp)    :: rij, coord_diff
 real(dp)     :: coord_sym1(3,3)
 real(dp)     :: coord_sym2(3,3)
 real(dp)     :: coord_temp(ndim,3)
 complex(dp) :: zi_vk_tn



zh(:,:) = cmplx(0.0_dp,0.0_dp,dp)

coord_temp(:,:)=coord(:,:)

coord_sym1(1,1)=(tn(1,2)+tn(1,3))/3.0_dp
coord_sym1(1,2)=(tn(2,2)+tn(2,3))/3.0_dp
coord_sym1(1,3)=coord_temp(ndim/2,3)

coord_sym1(2,1)=(tn(1,1)-tn(1,3))/3.0_dp
coord_sym1(2,2)=(tn(2,1)-tn(2,3))/3.0_dp
coord_sym1(2,3)=coord_temp(ndim/2,3)

coord_sym1(3,1)=-(tn(1,1)+tn(1,2))/3.0_dp
coord_sym1(3,2)=-(tn(2,1)+tn(2,2))/3.0_dp
coord_sym1(3,3)=coord_temp(ndim/2,3)

coord_sym2(1,1)=(tn(1,1)+tn(1,2))/3.0_dp
coord_sym2(1,2)=(tn(2,1)+tn(2,2))/3.0_dp
coord_sym2(1,3)=coord(ndim,3)

coord_sym2(2,1)=(tn(1,3)-tn(1,1))/3.0_dp
coord_sym2(2,2)=(tn(2,3)-tn(2,1))/3.0_dp
coord_sym2(2,3)=coord(ndim,3)

coord_sym2(3,1)=-(tn(1,3)+tn(1,2))/3.0_dp
coord_sym2(3,2)=-(tn(2,3)+tn(2,2))/3.0_dp
coord_sym2(3,3)=coord(ndim,3)


do l1=1,3
do l2=1,3

coord_temp(ndim/2,:)=coord_sym1(l1,:)
coord_temp(ndim,:)=coord_sym2(l2,:)


 do i = ndim/2,ndim,ndim/2
   do j = 1 , ndim

     rij = norm2(coord_temp(i,:)-coord_temp(j,:))

     coord_diff = (coord_temp(i,3)-coord_temp(j,3))**2

     if(coord_diff < deltaR .and. rij > deltaR) then
       zh(i,j) = zh(i,j)+ft(rij)
     endif
     if(coord_diff > deltaR)then
        zh(i,j) = zh(i,j)+ftperp(rij)
     endif


     do nt = 1 , 3
       rij = norm2(coord_temp(i,:)-[coord_temp(j,1:2)-tn(:,nt),coord_temp(j,3)])
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
       if(coord_diff < deltaR)then
         zh(i,j)=zh(i,j)+ft(rij)*exp( zi_vk_tn)
         zh(j,i)=zh(j,i)+ft(rij)*exp(-zi_vk_tn)
       else
         zh(i,j)=zh(i,j)+ftperp(rij)*exp( zi_vk_tn)
         zh(j,i)=zh(j,i)+ftperp(rij)*exp(-zi_vk_tn)
       endif
     end do
   end do
 end do

 do i = 1,ndim/2-1
    do j = ndim/2, ndim,ndim/2

      rij = norm2(coord_temp(i,:)-coord_temp(j,:))

      coord_diff = (coord_temp(i,3)-coord_temp(j,3))**2

      if(coord_diff < deltaR .and. rij > deltaR) then
        zh(i,j) = zh(i,j)+ft(rij)
      endif
       if(coord_diff > deltaR)then
         zh(i,j) = zh(i,j)+ftperp(rij)
      endif

 !!!!!! Coupling to adjacent Wigner-Seitz cell with tn
      do nt = 1 , 3
        rij = norm2(coord_temp(i,:)-[coord_temp(j,1:2)-tn(:,nt),coord_temp(j,3)])
        zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
        if(coord_diff < deltaR)then
          zh(i,j)=zh(i,j)+ft(rij)*exp( zi_vk_tn)
          zh(j,i)=zh(j,i)+ft(rij)*exp(-zi_vk_tn)
        else
          zh(i,j)=zh(i,j)+ftperp(rij)*exp( zi_vk_tn)
          zh(j,i)=zh(j,i)+ftperp(rij)*exp(-zi_vk_tn)
        endif
      end do
    end do
  end do

  do i = ndim/2+1,ndim-1
     do j = ndim/2, ndim,ndim/2

       rij = norm2(coord_temp(i,:)-coord_temp(j,:))

       coord_diff = (coord_temp(i,3)-coord_temp(j,3))**2

       if(coord_diff < deltaR .and. rij > deltaR) then
         zh(i,j) = zh(i,j)+ft(rij)
       endif
       if(coord_diff > deltaR)then
          zh(i,j) = zh(i,j)+ftperp(rij)
       endif

  !!!!!! Coupling to adjacent Wigner-Seitz cell with tn
       do nt = 1 , 3
         rij = norm2(coord_temp(i,:)-[coord_temp(j,1:2)-tn(:,nt),coord_temp(j,3)])
         zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
         if(coord_diff < deltaR)then
           zh(i,j)=zh(i,j)+ft(rij)*exp( zi_vk_tn)
           zh(j,i)=zh(j,i)+ft(rij)*exp(-zi_vk_tn)
         else
           zh(i,j)=zh(i,j)+ftperp(rij)*exp( zi_vk_tn)
           zh(j,i)=zh(j,i)+ftperp(rij)*exp(-zi_vk_tn)
         endif
       end do
     end do
   end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

enddo
enddo

zh(:,:)=zh(:,:)/cmplx(9.0_dp,0.0_dp,dp)

!!!!!!!!!!! Layer 1

 do i = 1 , ndim/2-1
   do j = 1 , ndim/2-1

     rij = norm2(coord_temp(i,:)-coord_temp(j,:))

       if(rij > deltaR) then
       zh(i,j) = zh(i,j)+ft(rij)
       endif

     do nt = 1 , 3
       rij = norm2(coord_temp(i,:)-[coord_temp(j,1:2)-tn(:,nt),coord_temp(j,3)])
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)

       zh(i,j)=zh(i,j)+ft(rij)*exp( zi_vk_tn)
       zh(j,i)=zh(j,i)+ft(rij)*exp(-zi_vk_tn)
     end do

   end do
 end do

!!!!!!!!!!! Layer 2

 do i = ndim/2+1,ndim-1
   do j = ndim/2+1,ndim-1

     rij = norm2(coord_temp(i,:)-coord_temp(j,:))

     if(rij > deltaR) then
       zh(i,j) = zh(i,j)+ft(rij)
     endif

     do nt = 1 , 3
       rij = norm2(coord_temp(i,:)-[coord_temp(j,1:2)-tn(:,nt),coord_temp(j,3)])
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)

       zh(i,j)=zh(i,j)+ft(rij)*exp( zi_vk_tn)
       zh(j,i)=zh(j,i)+ft(rij)*exp(-zi_vk_tn)
     end do

   end do
 end do

!!!!!!!!!!! Coupling Layer 1 to Layer 2

do i = 1 , ndim/2-1
   do j = ndim/2+1 , ndim-1

     rij = norm2(coord_temp(i,:)-coord_temp(j,:))

     zh(i,j) = zh(i,j)+ftperp(rij)

     do nt = 1 , 3
       rij = norm2(coord_temp(i,:)-[coord_temp(j,1:2)-tn(:,nt),coord_temp(j,3)])
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
         zh(i,j)=zh(i,j)+ftperp(rij)*exp( zi_vk_tn)
         zh(j,i)=zh(j,i)+ftperp(rij)*exp(-zi_vk_tn)
     end do

   end do
 end do

!!!!!!!!!!! Coupling Layer 2 to Layer 1

 do i = ndim/2+1 , ndim-1
    do j = 1 , ndim/2-1

      rij = norm2(coord_temp(i,:)-coord_temp(j,:))

      zh(i,j) = zh(i,j)+ftperp(rij)

      do nt = 1 , 3
        rij = norm2(coord_temp(i,:)-[coord_temp(j,1:2)-tn(:,nt),coord_temp(j,3)])
        zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
          zh(i,j)=zh(i,j)+ftperp(rij)*exp( zi_vk_tn)
          zh(j,i)=zh(j,i)+ftperp(rij)*exp(-zi_vk_tn)
      end do

    end do
  end do



end subroutine ComplexH0_sym


!c!c!cC
!c!c!cC
!c!c!cC

subroutine ComplexH0_Hartree_sym(zh,coord,potential,ndim,tn)

 integer(dp), intent(in)    :: ndim
 real(dp)   , intent(in)    :: tn(2,3), potential(ndim)
 real(dp)   , intent(in) :: coord(ndim,3)
 complex(dp), intent(inout) :: zh(ndim,ndim)

 integer(dp) :: i, j, nt, l1,l2
 real(dp)    :: rij, coord_diff
 real(dp)    :: coord_temp(ndim,3)
 real(dp)     :: coord_sym1(3,3)
 real(dp)     :: coord_sym2(3,3)


zh(:,:) = cmplx(0.0_dp,0.0_dp,dp)

coord_temp(:,:)=coord(:,:)

coord_sym1(1,1)=(tn(1,2)+tn(1,3))/3.0_dp
coord_sym1(1,2)=(tn(2,2)+tn(2,3))/3.0_dp
coord_sym1(1,3)=coord_temp(ndim/2,3)

coord_sym1(2,1)=(tn(1,1)-tn(1,3))/3.0_dp
coord_sym1(2,2)=(tn(2,1)-tn(2,3))/3.0_dp
coord_sym1(2,3)=coord_temp(ndim/2,3)

coord_sym1(3,1)=-(tn(1,1)+tn(1,2))/3.0_dp
coord_sym1(3,2)=-(tn(2,1)+tn(2,2))/3.0_dp
coord_sym1(3,3)=coord_temp(ndim/2,3)

coord_sym2(1,1)=(tn(1,1)+tn(1,2))/3.0_dp
coord_sym2(1,2)=(tn(2,1)+tn(2,2))/3.0_dp
coord_sym2(1,3)=coord_temp(ndim,3)

coord_sym2(2,1)=(tn(1,3)-tn(1,1))/3.0_dp
coord_sym2(2,2)=(tn(2,3)-tn(2,1))/3.0_dp
coord_sym2(2,3)=coord_temp(ndim,3)

coord_sym2(3,1)=-(tn(1,3)+tn(1,2))/3.0_dp
coord_sym2(3,2)=-(tn(2,3)+tn(2,2))/3.0_dp
coord_sym2(3,3)=coord_temp(ndim,3)


do l1=1,3
do l2=1,3

coord_temp(ndim/2,:)=coord_sym1(l1,:)
coord_temp(ndim,:)=coord_sym2(l2,:)


 do i = ndim/2,ndim,ndim/2
   do j = 1 , ndim

     rij = norm2(coord_temp(i,:)-coord_temp(j,:))

     coord_diff = (coord_temp(i,3)-coord_temp(j,3))**2

     if(coord_diff < deltaR) then
       zh(i,j) = zh(i,j)+ft(rij)
     else
        zh(i,j) = zh(i,j)+ftperp(rij)
     endif

   end do
 end do

 do i = 1,ndim/2-1
    do j = ndim/2, ndim,ndim/2

      rij = norm2(coord_temp(i,:)-coord_temp(j,:))

      coord_diff = (coord_temp(i,3)-coord_temp(j,3))**2

      if(coord_diff < deltaR) then
        zh(i,j) = zh(i,j)+ft(rij)
      else
         zh(i,j) = zh(i,j)+ftperp(rij)
      endif

    end do
  end do

  do i = ndim/2+1,ndim-1
     do j = ndim/2, ndim,ndim/2

       rij = norm2(coord_temp(i,:)-coord_temp(j,:))

       coord_diff = (coord_temp(i,3)-coord_temp(j,3))**2

       if(coord_diff < deltaR) then
         zh(i,j) = zh(i,j)+ft(rij)
       else
          zh(i,j) = zh(i,j)+ftperp(rij)
       endif

     end do
   end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



enddo
enddo

zh(:,:)=zh(:,:)/cmplx(9.0_dp,0.0_dp,dp)

!c!c!cC
!c!c!cC
!c!c!cC

do i = 1 , ndim
   zh(i,i)=zh(i,i)+potential(i)
enddo

!!!!!!!!!!! Layer 1

 do i = 1 , ndim/2-1
   do j = 1 , ndim/2-1

     rij = norm2(coord_temp(i,:)-coord_temp(j,:))

     zh(i,j) = zh(i,j)+ft(rij)

   end do
 end do

!!!!!!!!!!! Layer 2

 do i = ndim/2+1,ndim-1
   do j = ndim/2+1,ndim-1

     rij = norm2(coord_temp(i,:)-coord_temp(j,:))

     zh(i,j) = zh(i,j)+ft(rij)

   end do
 end do

!!!!!!!!!!! Coupling Layer 1 to Layer 2

do i = 1 , ndim/2-1
   do j = ndim/2+1 , ndim-1

     rij = norm2(coord_temp(i,:)-coord_temp(j,:))

     zh(i,j) = zh(i,j)+ftperp(rij)

   end do
 end do

!!!!!!!!!!! Coupling Layer 2 to Layer 1

 do i = ndim/2+1 , ndim-1
    do j = 1 , ndim/2-1

      rij = norm2(coord_temp(i,:)-coord_temp(j,:))

      zh(i,j) = zh(i,j)+ftperp(rij)

    end do
  end do


end subroutine ComplexH0_Hartree_sym

!c!c!cC
!c!c!cC
!c!c!cC

subroutine ComplexHk_Hartree_sym(zh,coord,ndim,vk,tn)

 integer(dp), intent(in)    :: ndim
 real(dp)   , intent(in)    :: tn(2,3), vk(2)
 real(dp)   , intent(in) :: coord(ndim,3)
 complex(dp), intent(inout) :: zh(ndim,ndim)

 integer(dp) :: i, j, nt, l1,l2
 real(dp)    :: rij, coord_diff
 real(dp)     :: coord_sym1(3,3)
 real(dp)     :: coord_sym2(3,3)
 real(dp)     :: coord_temp(ndim,3)
 complex(dp) :: zi_vk_tn


zh(:,:) = cmplx(0.0_dp,0.0_dp,dp)

coord_temp(:,:)=coord(:,:)

coord_sym1(1,1)=(tn(1,2)+tn(1,3))/3.0_dp
coord_sym1(1,2)=(tn(2,2)+tn(2,3))/3.0_dp
coord_sym1(1,3)=coord_temp(ndim/2,3)

coord_sym1(2,1)=(tn(1,1)-tn(1,3))/3.0_dp
coord_sym1(2,2)=(tn(2,1)-tn(2,3))/3.0_dp
coord_sym1(2,3)=coord_temp(ndim/2,3)

coord_sym1(3,1)=-(tn(1,1)+tn(1,2))/3.0_dp
coord_sym1(3,2)=-(tn(2,1)+tn(2,2))/3.0_dp
coord_sym1(3,3)=coord_temp(ndim/2,3)

coord_sym2(1,1)=(tn(1,1)+tn(1,2))/3.0_dp
coord_sym2(1,2)=(tn(2,1)+tn(2,2))/3.0_dp
coord_sym2(1,3)=coord_temp(ndim,3)

coord_sym2(2,1)=(tn(1,3)-tn(1,1))/3.0_dp
coord_sym2(2,2)=(tn(2,3)-tn(2,1))/3.0_dp
coord_sym2(2,3)=coord_temp(ndim,3)

coord_sym2(3,1)=-(tn(1,3)+tn(1,2))/3.0_dp
coord_sym2(3,2)=-(tn(2,3)+tn(2,2))/3.0_dp
coord_sym2(3,3)=coord_temp(ndim,3)


do l1=1,3
do l2=1,3

coord_temp(ndim/2,:)=coord_sym1(l1,:)
coord_temp(ndim,:)=coord_sym2(l2,:)


 do i = ndim/2,ndim,ndim/2
   do j = 1 , ndim

     coord_diff = (coord_temp(i,3)-coord_temp(j,3))**2

     do nt = 1 , 3
       rij = norm2(coord_temp(i,:)-[coord_temp(j,1:2)-tn(:,nt),coord_temp(j,3)])
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
       if(coord_diff < deltaR)then
         zh(i,j)=zh(i,j)+ft(rij)*exp( zi_vk_tn)
         zh(j,i)=zh(j,i)+ft(rij)*exp(-zi_vk_tn)
       else
         zh(i,j)=zh(i,j)+ftperp(rij)*exp( zi_vk_tn)
         zh(j,i)=zh(j,i)+ftperp(rij)*exp(-zi_vk_tn)
       endif
     end do

   end do
 end do

 do i = 1,ndim/2-1
    do j = ndim/2, ndim,ndim/2

      coord_diff = (coord_temp(i,3)-coord_temp(j,3))**2

 !!!!!! Coupling to adjacent Wigner-Seitz cell with tn
      do nt = 1 , 3
        rij = norm2(coord_temp(i,:)-[coord_temp(j,1:2)-tn(:,nt),coord_temp(j,3)])
        zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
        if(coord_diff < deltaR)then
          zh(i,j)=zh(i,j)+ft(rij)*exp( zi_vk_tn)
          zh(j,i)=zh(j,i)+ft(rij)*exp(-zi_vk_tn)
        else
          zh(i,j)=zh(i,j)+ftperp(rij)*exp( zi_vk_tn)
          zh(j,i)=zh(j,i)+ftperp(rij)*exp(-zi_vk_tn)
        endif
      end do

    end do
  end do

  do i = ndim/2+1,ndim-1
     do j = ndim/2, ndim,ndim/2

       coord_diff = (coord_temp(i,3)-coord_temp(j,3))**2

  !!!!!! Coupling to adjacent Wigner-Seitz cell with tn
       do nt = 1 , 3
         rij = norm2(coord_temp(i,:)-[coord_temp(j,1:2)-tn(:,nt),coord_temp(j,3)])
         zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
         if(coord_diff < deltaR)then
           zh(i,j)=zh(i,j)+ft(rij)*exp( zi_vk_tn)
           zh(j,i)=zh(j,i)+ft(rij)*exp(-zi_vk_tn)
         else
           zh(i,j)=zh(i,j)+ftperp(rij)*exp( zi_vk_tn)
           zh(j,i)=zh(j,i)+ftperp(rij)*exp(-zi_vk_tn)
         endif
       end do

     end do
   end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



enddo
enddo

zh(:,:)=zh(:,:)/cmplx(9.0_dp,0.0_dp,dp)

!!!!!!!!!!! Layer 1

 do i = 1 , ndim/2-1
   do j = 1 , ndim/2-1

     do nt = 1 , 3
       rij = norm2(coord_temp(i,:)-[coord_temp(j,1:2)-tn(:,nt),coord_temp(j,3)])
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)

       zh(i,j)=zh(i,j)+ft(rij)*exp( zi_vk_tn)
       zh(j,i)=zh(j,i)+ft(rij)*exp(-zi_vk_tn)
     end do

   end do
 end do

!!!!!!!!!!! Layer 2

 do i = ndim/2+1,ndim-1
   do j = ndim/2+1,ndim-1

     do nt = 1 , 3
       rij = norm2(coord_temp(i,:)-[coord_temp(j,1:2)-tn(:,nt),coord_temp(j,3)])
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)

       zh(i,j)=zh(i,j)+ft(rij)*exp( zi_vk_tn)
       zh(j,i)=zh(j,i)+ft(rij)*exp(-zi_vk_tn)
     end do

   end do
 end do

!!!!!!!!!!! Coupling Layer 1 to Layer 2

do i = 1 , ndim/2-1
   do j = ndim/2+1 , ndim-1

     do nt = 1 , 3
       rij = norm2(coord_temp(i,:)-[coord_temp(j,1:2)-tn(:,nt),coord_temp(j,3)])
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
         zh(i,j)=zh(i,j)+ftperp(rij)*exp( zi_vk_tn)
         zh(j,i)=zh(j,i)+ftperp(rij)*exp(-zi_vk_tn)
     end do

   end do
 end do

!!!!!!!!!!! Coupling Layer 2 to Layer 1

 do i = ndim/2+1 , ndim-1
    do j = 1 , ndim/2-1

      do nt = 1 , 3
        rij = norm2(coord_temp(i,:)-[coord_temp(j,1:2)-tn(:,nt),coord_temp(j,3)])
        zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
          zh(i,j)=zh(i,j)+ftperp(rij)*exp( zi_vk_tn)
          zh(j,i)=zh(j,i)+ftperp(rij)*exp(-zi_vk_tn)
      end do

    end do
  end do


end subroutine ComplexHk_Hartree_sym


!c!c!cC
!c!c!cC
!c!c!cC




!c!c!cC
!c!c!cC
!c!c!cC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure function fv_e(r)
real(dp)   , intent(in)    :: r
real(dp) :: fv_e

real(dp) x,I0,K0
real(dp) xi
integer(dp) n

xi=100.0_dp/2.46_dp/3._dp ! 10nm/a

fv_e=0._dp

if(r.NE.0)then

if(r.LT.xi/3.)then

do n=0,floor((2._dp-pi*(r/xi))/2._dp/pi/(r/xi))

x=(2*n+1)*pi*(r/xi)

I0=cosh(x)/(1 + 0.5**2*x**2._dp)**0.25*(1 + 0.24273*x**2._dp)/(1 + 0.43023*x**2._dp)

K0=-log(x/2._dp)*I0-0.5772156648942439+0.42278433434244916*(x/2._dp)**2+0.23069609660563425*(x/2._dp)**4 &
+0.03489207637875737*(x/2._dp)**6+0.002615030023757213*(x/2._dp)**8+0.00011811080908871537*(x/2._dp)**10+0.000003889449474816304*(x/2._dp)**12

fv_e=fv_e+4.*K0

enddo

do n=floor((2._dp-pi*(r/xi))/2._dp/pi/(r/xi))+1,100

x=(2*n+1)*pi*(r/xi)

K0=(1.2533127470318168-0.07830516193156768*(2._dp/x)+0.021807436132174653*(2._dp/x)**2 &
-0.010428688609726261*(2._dp/x)**3 + 0.005672414173632901*(2._dp/x)**4-0.0024265259192435785*(2._dp/x)**5 &
+0.0005249625381161658*(2._dp/x)**6)*exp(-x)/sqrt(x)

fv_e=fv_e+4.*K0

enddo

fv_e=fv_e/xi

else

fv_e=2.0_dp*sqrt(2.0_dp)/sqrt(r/xi)*exp(-pi*r/xi)/xi

endif

else
fv_e=0.0_dp
endif

return
end function fv_e

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




pure function fv_a(r)
real(dp)   , intent(in)    :: r
real(dp) :: fv_a

real(dp) xi,alpha0,t0,epsilon0

alpha0=2.5_dp
t0=2.7_dp
epsilon0=6.0_dp
xi=100.0_dp/2.46_dp/3._dp ! 10nm/a

fv_a=2.0_dp*sqrt(2.0_dp)/sqrt(r/xi)*exp(-pi*r/xi)/xi

if(r.EQ.0)then
fv_a=0.0_dp
endif

return
end function fv_a


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure function ftperp(r)

 real(dp), intent(in) :: r
 real(dp) :: ftperp

 real(dp) :: pd

 pd=(d0/r)**2

 ftperp = (-2.7_dp*exp((a0-r)/r0)*(1.0_dp-pd)+0.48_dp*exp((d0*pressure-r)/r0)*pd)*tperp

end function ftperp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure function ft(r)

 real(dp), intent(in) :: r
 real(dp) :: ft

 !if( r <= 1.01_dp*a0.AND.r>deltaR) then
 if( r>deltaR) then
   ft = -2.7_dp*exp((a0-r)/r0)
 else
   ft = 0.0_dp
 end if

end function ft

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine vx_new(zvx,ndim,coord,vk,tn,xy)

 integer(dp), intent(in)    :: ndim, xy
 real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3)
 complex(dp), intent(inout) :: zvx(:,:)

 integer(dp) :: i, j, nt
 real(dp) :: rij, rxy, coord_diff
 complex(dp) :: zi_vk_tn

 zvx = cmplx(0.0_dp,0.0_dp,dp)

 do i = 1 , ndim
   do j = i , ndim

     rij = norm2(coord(i,:)-coord(j,:))
     rxy  = coord(i,xy)-coord(j,xy)
     coord_diff = (coord(i,3)-coord(j,3))**2

     if(coord_diff < deltaR .and. rij > deltaR) then
       zvx(i,j) = cmplx(0.0_dp,ft(rij)*rxy,dp)
     end if

     if(coord_diff > deltaR)then
       zvx(i,j) = cmplx(0.0_dp,ftperp(rij)*rxy,dp)
     endif

!!!!!!! Coupling to adjacent Wigner-Seitz cell with tn
     !do nt = 1 , 3
      ! rij = norm2(coord(i,:)-[coord(j,1:2)-tn(:,nt),coord(j,3)])
       !rxy=coord(i,xy)-coord(j,xy)+tn(xy,nt)
       !zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)

       !if(coord_diff < deltaR) then
        ! zvx(i,j)=zvx(i,j)+ft(rij)*exp( zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
        ! zvx(j,i)=zvx(j,i)-ft(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
       !else
        ! zvx(i,j)=zvx(i,j)+ftperp(rij)*exp( zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
        ! zvx(j,i)=zvx(j,i)-ftperp(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
       !endif
     !end do

     !!!!!! Coupling to adjacent Wigner-Seitz cell with tn
     do nt = 1 , 3
       rij = norm2(coord(i,:)-[coord(j,1:2)-tn(:,nt),coord(j,3)])
       rxy=coord(i,xy)-coord(j,xy)+tn(xy,nt)
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
       if(coord_diff < deltaR)then
         zvx(i,j)=zvx(i,j)+ft(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
       else
         zvx(i,j)=zvx(i,j)+ftperp(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
       endif
     !!!!!! Coupling to adjacent Wigner-Seitz cell with -tn
    rij = norm2(coord(i,:)-[coord(j,1:2)+tn(:,nt),coord(j,3)])
    rxy=coord(i,xy)-coord(j,xy)-tn(xy,nt)
    zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,tn(:,nt)),dp)
    if(coord_diff < deltaR)then
        zvx(i,j)=zvx(i,j)+ft(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
    else
        zvx(i,j)=zvx(i,j)+ftperp(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
    endif
    enddo

   end do
 end do

end subroutine vx_new

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine vxy12I(zvxy1,zvxy2,zvxyI,ndim,coord,vk,tn,xy)

 integer(dp), intent(in)    :: ndim, xy
 real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3)
 complex(dp), intent(inout) :: zvxy1(:,:),zvxy2(:,:),zvxyI(:,:)

 integer(dp) :: i, j, nt
 real(dp) :: rij, rxy, coord_diff
 complex(dp) :: zi_vk_tn

 zvxy1 = cmplx(0.0_dp,0.0_dp,dp)
 zvxy2 = cmplx(0.0_dp,0.0_dp,dp)
 zvxyI = cmplx(0.0_dp,0.0_dp,dp)

 do i = 1 , ndim
   do j = i , ndim

     rij = norm2(coord(i,:)-coord(j,:))
     rxy  = coord(i,xy)-coord(j,xy)

     if(     coord(i,3)<deltaR .and. coord(j,3)<deltaR .and. rij>deltaR) then
       zvxy1(i,j) = cmplx(0.0_dp,ft(rij)*rxy,dp)
     else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR .and. rij>deltaR) then
       zvxy2(i,j) = cmplx(0.0_dp,ft(rij)*rxy,dp)
     else
       zvxyI(i,j) = cmplx(0.0_dp,ftperp(rij)*rxy,dp)
     endif

     !!!!!! Coupling to adjacent Wigner-Seitz cell with tn
     do nt = 1 , 3
       rij = norm2(coord(i,:)-[coord(j,1:2)-tn(:,nt),coord(j,3)])
       rxy=coord(i,xy)-coord(j,xy)+tn(xy,nt)
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
       if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
         zvxy1(i,j)=zvxy1(i,j)+ft(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
       else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
         zvxy2(i,j)=zvxy2(i,j)+ft(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
       else
         zvxyI(i,j)=zvxyI(i,j)+ftperp(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
       endif
     !!!!!! Coupling to adjacent Wigner-Seitz cell with -tn
    rij = norm2(coord(i,:)-[coord(j,1:2)+tn(:,nt),coord(j,3)])
    rxy=coord(i,xy)-coord(j,xy)-tn(xy,nt)
    zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,tn(:,nt)),dp)
    if(coord(i,3)<deltaR .and. coord(j,3)<deltaR)then
        zvxy1(i,j)=zvxy1(i,j)+ft(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
    else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
        zvxy2(i,j)=zvxy2(i,j)+ft(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
        else
         zvxyI(i,j)=zvxyI(i,j)+ftperp(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
    endif
    enddo

   end do
 end do

end subroutine vxy12I

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine vxy1_new(zvxy1,zvxy2,ndim,coord,vk,tn,xy)

 integer(dp), intent(in)    :: ndim, xy
 real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3)
 complex(dp), intent(inout) :: zvxy1(:,:), zvxy2(:,:)

 integer(dp) :: i, j, nt
 real(dp) :: rij, rxy
 complex(dp) :: zi_vk_tn

 zvxy1 = cmplx(0.0_dp,0.0_dp,dp)
 zvxy2 = cmplx(0.0_dp,0.0_dp,dp)

 do i = 1 , ndim
   do j = 1 , ndim

     rij = norm2(coord(i,:)-coord(j,:))
     rxy = coord(i,xy)-coord(j,xy)

     if(     coord(i,3)<deltaR .and. coord(j,3)<deltaR .and. rij>deltaR) then
       zvxy1(i,j) = cmplx(0.0_dp,ft(rij)*rxy,dp)
     else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR .and. rij>deltaR) then
       zvxy2(i,j) = cmplx(0.0_dp,ft(rij)*rxy,dp)
     endif

!!!!!!! Coupling to adjacent Wigner-Seitz cell with tn
     do nt = 1 , 3
       rij = norm2(coord(i,:)-[coord(j,1:2)-tn(:,nt),coord(j,3)])
       rxy = (coord(i,xy)-coord(j,xy)+tn(xy,nt))
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)

       if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
         zvxy1(i,j)=zvxy1(i,j)+ft(rij)*exp( zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
         zvxy1(j,i)=zvxy1(j,i)-ft(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
       else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
         zvxy2(i,j)=zvxy2(i,j)+ft(rij)*exp( zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
         zvxy2(j,i)=zvxy2(j,i)-ft(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
       endif
     end do

   end do
 end do

end subroutine vxy1_new

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine vxyI_new(zvxy,ndim,coord,vk,tn,xy)

 integer(dp), intent(in)    :: ndim, xy
 real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3)
 complex(dp), intent(inout) :: zvxy(:,:)

 integer(dp) :: i, j, nt
 real(dp) :: rij, rxy, coord_diff
 complex(dp) :: zi_vk_tn

 zvxy = cmplx(0.0_dp,0.0_dp,dp)

 do i = 1 , ndim
   do j = 1 , ndim

     rij = norm2(coord(i,:)-coord(j,:))
     rxy = coord(i,xy)-coord(j,xy)
     coord_diff = (coord(i,3)-coord(j,3))**2

     if(coord_diff>deltaR) then
       zvxy(i,j) = cmplx(0.0_dp,ftperp(rij)*rxy,dp)
     endif

!!!!!!! Coupling to adjacent Wigner-Seitz cell with tn
     do nt = 1 , 3
       rij = norm2(coord(i,:)-[coord(j,1:2)-tn(:,nt),coord(j,3)])
       rxy=(coord(i,xy)-coord(j,xy)+tn(1, nt))
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)

       if(coord_diff>deltaR) then
         zvxy(i,j)=zvxy(i,j)+ftperp(rij)*exp( zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
         zvxy(j,i)=zvxy(j,i)-ftperp(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)
       endif
     end do

   end do
 end do

end subroutine vxyI_new

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine vxy1D_new(zvxy,ndim,coord,vk,tn,xy)

 integer(dp), intent(in)    :: ndim, xy
 real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3)
 complex(dp), intent(inout) :: zvxy(:,:)

 integer(dp) :: i, j, nt
 real(dp) :: rij, rxy, coord_diff
 complex(dp) :: zi_vk_tn

 zvxy = cmplx(0.0_dp,0.0_dp,dp)

 do i = 1 , ndim
   do j = 1 , ndim

     rij = norm2(coord(i,:)-coord(j,:))
     rxy = coord(i,xy)-coord(j,xy)
     coord_diff = (coord(i,3)-coord(j,3))**2

     if(coord(i,3)<deltaR .and. coord(j,3)<deltaR .and. rij>deltaR) then
       zvxy(i,j)=ft(rij)*rxy**2
     endif

!!!!!!! Coupling to adjacent Wigner-Seitz cell with tn
     do nt = 1 , 3
       rij = norm2(coord(i,:)-[coord(j,1:2)-tn(:,nt),coord(j,3)])
       rxy=(coord(i,xy)-coord(j,xy)+tn(1, nt))
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)

       if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
         zvxy(i,j)=zvxy(i,j)+ft(rij)*exp( zi_vk_tn)*rxy**2
         zvxy(j,i)=zvxy(j,i)+ft(rij)*exp(-zi_vk_tn)*rxy**2
       endif
     end do

   end do
 end do


end subroutine vxy1D_new

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!c!c!cC
!c!c!cC
!c!c!cC

subroutine vxy12ID_MF(zvxy1,zvxy2,zvxyI,coord,alpha,zfock,ind_j,ind_l,ndim,ind,vk,tn,nxy)

integer(dp), intent(in)    :: ndim, ind, nxy

integer(dp), intent(in)    :: ind_j(ind), ind_l(ind)
real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3),alpha
complex(dp) , intent(in)   :: zfock(ndim,ndim,ind)
complex(dp), intent(out) :: zvxy1(:,:),zvxy2(:,:),zvxyI(:,:)

integer(dp) :: i, j, nt, icount, n1, n2
real(dp)    :: rij, rxy, coord_diff
complex(dp) :: zi_vk_tn,zsum1,zsum2,zsumI,zphase


 zvxy1 = cmplx(0.0_dp,0.0_dp,dp)
 zvxy2 = cmplx(0.0_dp,0.0_dp,dp)
 zvxyI = cmplx(0.0_dp,0.0_dp,dp)


 do i = 1 , ndim
   do j = i , ndim

     zphase=exp(cmplx(0.0_dp,-dot_product(vk,coord(i,1:2)-coord(j,1:2))*fphase,dp))

     rij = norm2(coord(i,:)-coord(j,:))
     rxy = coord(i,nxy)-coord(j,nxy)


if(coord(i,3)<deltaR .and. coord(j,3)<deltaR .and. rij>deltaR) then
  zvxy1(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*rxy*rxy*0.5_dp
  zvxy1(i,j)=zvxy1(i,j)+ft(rij)*zphase*rxy*rxy
else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR .and. rij>deltaR) then
  zvxy2(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*rxy*rxy*0.5_dp
  zvxy2(i,j)=zvxy2(i,j)+ft(rij)*zphase*rxy*rxy
else if(coord(i,3).NE.coord(j,3)) then
  zvxyI(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*rxy*rxy*0.5_dp
  zvxyI(i,j)=zvxyI(i,j) + ftperp(rij)*zphase*rxy*rxy
endif

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

     do icount=2,ind
     n1=ind_j(icount)
     n2=ind_l(icount)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
      zvxy1(i,j)=zvxy1(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*rxy*rxy
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
      zvxy2(i,j)=zvxy2(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*rxy*rxy
   else if(coord(i,3).NE.coord(j,3)) then
      zvxyI(i,j)=zvxyI(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase*rxy*rxy
   endif

   rij = norm2(coord(i,:)-[coord(j,1:2)+n1*tn(:,1)+n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)-n1*tn(nxy,1)-n2*tn(nxy,2)

   zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
      zvxy1(i,j)=zvxy1(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*rxy*rxy
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
      zvxy2(i,j)=zvxy2(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*rxy*rxy
   else if(coord(i,3).NE.coord(j,3)) then
      zvxyI(i,j)=zvxyI(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase*rxy*rxy
   endif

   enddo

!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!

   zsum1=cmplx(0.0_dp,0.0_dp,dp)
   zsum2=cmplx(0.0_dp,0.0_dp,dp)
   zsumI=cmplx(0.0_dp,0.0_dp,dp)
   do icount=2,ind
   n1=ind_j(icount)
   n2=ind_l(icount)

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
      zsum1=zsum1-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*rxy*rxy
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
      zsum2=zsum2-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*rxy*rxy
   else if(coord(i,3).NE.coord(j,3)) then
      zsumI=zsumI-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*rxy*rxy
   endif

   rij = norm2(coord(j,:)-[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])
   rxy=coord(j,nxy)-coord(i,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
   zsum1=zsum1+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*rxy*rxy
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
   zsum2=zsum2+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*rxy*rxy
   else if(coord(i,3).NE.coord(j,3)) then
   zsumI=zsumI+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*rxy*rxy
   endif

   enddo

   zvxy1(i,j)=zvxy1(i,j)+zsum1*zphase*0.5_dp
   zvxy2(i,j)=zvxy2(i,j)+zsum2*zphase*0.5_dp
   zvxyI(i,j)=zvxyI(i,j)+zsumI*zphase*0.5_dp

   end do
 end do



end subroutine vxy12ID_MF

!c!c!cC
!c!c!cC
!c!c!cC

subroutine vxy12ID_NO(zvxy1,zvxy2,zvxyI,coord,alpha,zfock,ind_j,ind_l,ndim,ind,vk,tn,nxy)

integer(dp), intent(in)    :: ndim, ind, nxy

integer(dp), intent(in)    :: ind_j(ind), ind_l(ind)
real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3),alpha
complex(dp) , intent(in)   :: zfock(ndim,ndim,ind)
complex(dp), intent(out) :: zvxy1(:,:),zvxy2(:,:),zvxyI(:,:)

integer(dp) :: i, j, nt, icount, n1, n2
real(dp)    :: rij, rxy, coord_diff
complex(dp) :: zi_vk_tn,zsum1,zsum2,zsumI,zphase


 zvxy1 = cmplx(0.0_dp,0.0_dp,dp)
 zvxy2 = cmplx(0.0_dp,0.0_dp,dp)
 zvxyI = cmplx(0.0_dp,0.0_dp,dp)


 do i = 1 , ndim
   do j = i , ndim

     zphase=exp(cmplx(0.0_dp,-dot_product(vk,coord(i,1:2)-coord(j,1:2))*fphase,dp))

     rij = norm2(coord(i,:)-coord(j,:))
     rxy = coord(i,nxy)-coord(j,nxy)


if(coord(i,3)<deltaR .and. coord(j,3)<deltaR .and. rij>deltaR) then
  zvxy1(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*rxy*rxy
else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR .and. rij>deltaR) then
  zvxy2(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*rxy*rxy
else if(coord(i,3).NE.coord(j,3)) then
  zvxyI(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*rxy*rxy
endif


!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!

   zsum1=cmplx(0.0_dp,0.0_dp,dp)
   zsum2=cmplx(0.0_dp,0.0_dp,dp)
   zsumI=cmplx(0.0_dp,0.0_dp,dp)
   do icount=2,ind
   n1=ind_j(icount)
   n2=ind_l(icount)

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
      zsum1=zsum1-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*rxy*rxy
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
      zsum2=zsum2-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*rxy*rxy
   else if(coord(i,3).NE.coord(j,3)) then
      zsumI=zsumI-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*rxy*rxy
   endif

   rij = norm2(coord(j,:)-[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])
   rxy=coord(j,nxy)-coord(i,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
   zsum1=zsum1+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*rxy*rxy
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
   zsum2=zsum2+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*rxy*rxy
   else if(coord(i,3).NE.coord(j,3)) then
   zsumI=zsumI+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*rxy*rxy
   endif

   enddo

   zvxy1(i,j)=zvxy1(i,j)+zsum1*zphase
   zvxy2(i,j)=zvxy2(i,j)+zsum2*zphase
   zvxyI(i,j)=zvxyI(i,j)+zsumI*zphase

   end do
 end do



end subroutine vxy12ID_NO

!c!c!cC
!c!c!cC
!c!c!cC

subroutine vxy12D_MF(zvxy,coord,alpha,zfock,ind_j,ind_l,ndim,ind,vk,tn,nxy)

integer(dp), intent(in)    :: ndim, ind, nxy

integer(dp), intent(in)    :: ind_j(ind), ind_l(ind)
real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3),alpha
complex(dp) , intent(in)   :: zfock(ndim,ndim,ind)
complex(dp), intent(out) :: zvxy(:,:)

integer(dp) :: i, j, nt, icount, n1, n2
real(dp)    :: rij, rxy, coord_diff
complex(dp) :: zi_vk_tn,zsum1,zsum2,zsumI,zphase


 zvxy = cmplx(0.0_dp,0.0_dp,dp)


 do i = 1 , ndim
   do j = i , ndim

     zphase=exp(cmplx(0.0_dp,-dot_product(vk,coord(i,1:2)-coord(j,1:2))*fphase,dp))

     rij = norm2(coord(i,:)-coord(j,:))
     rxy = coord(i,nxy)-coord(j,nxy)


if(coord(i,3)<deltaR .and. coord(j,3)<deltaR .and. rij>deltaR) then
  zvxy(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*rxy*rxy
  zvxy(i,j)=zvxy(i,j)+ft(rij)*zphase*rxy*rxy
else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR .and. rij>deltaR) then
  zvxy(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*rxy*rxy
  zvxy(i,j)=zvxy(i,j)+ft(rij)*zphase*rxy*rxy
else if(coord(i,3).NE.coord(j,3)) then
  zvxy(i,j)=-alpha*fv(coord(i,:),coord(j,:))*conjg(zfock(i,j,1))*zphase*rxy*rxy
  zvxy(i,j)=zvxy(i,j) + ftperp(rij)*zphase*rxy*rxy
endif

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

     do icount=2,ind
     n1=ind_j(icount)
     n2=ind_l(icount)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
      zvxy(i,j)=zvxy(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*rxy*rxy
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
      zvxy(i,j)=zvxy(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*rxy*rxy
   else if(coord(i,3).NE.coord(j,3)) then
      zvxy(i,j)=zvxy(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase*rxy*rxy
   endif

   rij = norm2(coord(i,:)-[coord(j,1:2)+n1*tn(:,1)+n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)-n1*tn(nxy,1)-n2*tn(nxy,2)

   zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
      zvxy(i,j)=zvxy(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*rxy*rxy
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
      zvxy(i,j)=zvxy(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase*rxy*rxy
   else if(coord(i,3).NE.coord(j,3)) then
      zvxy(i,j)=zvxy(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase*rxy*rxy
   endif

   enddo

!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!

   zsum1=cmplx(0.0_dp,0.0_dp,dp)
   zsum2=cmplx(0.0_dp,0.0_dp,dp)
   zsumI=cmplx(0.0_dp,0.0_dp,dp)
   do icount=2,ind
   n1=ind_j(icount)
   n2=ind_l(icount)

   zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)

   rij = norm2(coord(i,:)-[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])
   rxy=coord(i,nxy)-coord(j,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
      zsum1=zsum1-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*rxy*rxy
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
      zsum2=zsum2-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*rxy*rxy
   else if(coord(i,3).NE.coord(j,3)) then
      zsumI=zsumI-alpha*fv(coord(i,:),[coord(j,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(j,3)])*conjg(zfock(i,j,icount))*exp(-zi_vk_tn)*rxy*rxy
   endif

   rij = norm2(coord(j,:)-[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])
   rxy=coord(j,nxy)-coord(i,nxy)+n1*tn(nxy,1)+n2*tn(nxy,2)

   if(coord(i,3)<deltaR .and. coord(j,3)<deltaR) then
   zsum1=zsum1+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*rxy*rxy
   else if(coord(i,3)>deltaR .and. coord(j,3)>deltaR) then
   zsum2=zsum2+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*rxy*rxy
   else if(coord(i,3).NE.coord(j,3)) then
   zsumI=zsumI+alpha*fv(coord(j,:),[coord(i,1:2)-n1*tn(:,1)-n2*tn(:,2),coord(i,3)])*zfock(j,i,icount)*exp( zi_vk_tn)*rxy*rxy
   endif

   enddo

   zvxy(i,j)=zvxy(i,j)+zsum1*zphase+zsum2*zphase+zsumI*zphase

   end do
 end do



end subroutine vxy12D_MF


subroutine vxyID_new(zvxy,ndim,coord,vk,tn,xy)

 integer(dp), intent(in)    :: ndim, xy
 real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3)
 complex(dp), intent(inout) :: zvxy(:,:)

 integer(dp) :: i, j, nt
 real(dp) :: rij, rxy, coord_diff
 complex(dp) :: zi_vk_tn

 zvxy = cmplx(0.0_dp,0.0_dp,dp)

 do i = 1 , ndim
   do j = 1 , ndim

     rij = norm2(coord(i,:)-coord(j,:))
     rxy = coord(i,xy)-coord(j,xy)
     coord_diff = (coord(i,3)-coord(j,3))**2

     if(coord_diff>deltaR) then
       zvxy(i,j)=ftperp(rij)*rxy**2
     endif

!!!!!!! Coupling to adjacent Wigner-Seitz cell with tn
     do nt = 1 , 3
       rij = norm2(coord(i,:)-[coord(j,1:2)-tn(:,nt),coord(j,3)])
       rxy=(coord(i,xy)-coord(j,xy)+tn(1, nt))
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)

       if(coord_diff>deltaR) then
         zvxy(i,j)=zvxy(i,j)+ftperp(rij)*exp( zi_vk_tn)*rxy**2
         zvxy(j,i)=zvxy(j,i)+ftperp(rij)*exp(-zi_vk_tn)*rxy**2
       endif
     end do

   end do
 end do

end subroutine vxyID_new

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
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

subroutine plotBands(bandsTB,fmu,npointsBack,ndimEV,numk,ncount,vq1,vq12)

 integer(dp), intent(in)    :: ndimEV,numk,ncount
 integer(dp), intent(in)    :: npointsBack(numk,numk)
 real(dp)   , intent(in)    :: fmu,vq1(2), vq12(2), bandsTB(ndimEV,ncount)

 integer(dp) :: i,ivk,ivk1,ivk2
 real(dp) :: vK1x,vK1y,vK2x,vK2y,vM1x,vM1y,vM2x,vM2y
 real(dp) :: a11,a12,a21,a22,b11,b12,b21,b22,vkx,vky,det
 real(dp) :: aGM, aKG, aMK, aT

 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 vK1x=(vq1(1)+vq12(1))/3.0_dp
 vK1y=(vq1(2)+vq12(2))/3.0_dp

 vK2x=2.0_dp*(vq1(1)+vq12(1))/3.0_dp
 vK2y=2.0_dp*(vq1(2)+vq12(2))/3.0_dp

 vM1x=vq1(1)/2.0_dp
 vM1y=vq1(2)/2.0_dp

 vM2x=vq12(1)/2.0_dp
 vM2y=vq12(2)/2.0_dp


 aKG=sqrt(vK1x**2+vK1y**2)
 aGM=sqrt(vM1x**2+vM1y**2)
 aMK=sqrt((vK1x-vM1x)**2+(vK1y-vM1y)**2)
 aT=aKG+aGM+aMK

 a11=vq1(1)/real(numk,dp)
 a12=vq12(1)/real(numk,dp)
 a21=vq1(2)/real(numk,dp)
 a22=vq12(2)/real(numk,dp)

 det=a11*a22-a12*a21
 b11=a22/det
 b12=-a12/det
 b21=-a21/det
 b22=a11/det

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

 do ivk2=numk/3,0,-1
     ivk1=ivk2
     !write(*,*) (numk/3-ivk2)*aKG/real(numk/3,dp)/aT,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
     write(99,fmt='(F16.8,3X)', advance="no") (numk/3-ivk2)*aKG/real(numk/3,dp)/aT
     do i=1,ndimEV-1
     write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(ivk1+1,ivk2+1))-fmu
     enddo
     write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(ivk1+1,ivk2+1))-fmu
 enddo
 do ivk2=1,numk/2
     ivk1=0
     !write(*,*) (aKG+(ivk2)*aGM/real(numk/2,dp))/aT,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
     write(99,fmt='(F16.8,3X)', advance="no") (aKG+(ivk2)*aGM/real(numk/2,dp))/aT
     do i=1,ndimEV-1
     write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(ivk1+1,ivk2+1))-fmu
     enddo
     write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(ivk1+1,ivk2+1))-fmu
 enddo
 do ivk=1,numk/6
     vkx=vM1x+(vK1x-vM1x)*ivk/real(numk/6,dp)
     vky=vM1y+(vK1y-vM1y)*ivk/real(numk/6,dp)
     ivk2=nint(b11*vkx+b12*vky)
     ivk1=nint(b21*vkx+b22*vky)
     !write(*,*) (aKG+aGM+(ivk)*aMK/real(numk/6,dp))/aT,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
     write(99,fmt='(F16.8,3X)', advance="no") (aKG+aGM+(ivk)*aMK/real(numk/6,dp))/aT
     do i=1,ndimEV-1
     write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(ivk1+1,ivk2+1))-fmu
     enddo
     write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(ivk1+1,ivk2+1))-fmu
 enddo

 close(99)

end subroutine plotBands

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine plotBandsGamma(bandsTB,fmu,npointsBack,ndimEV,numk,ncount,vq1,vq12)

 integer(dp), intent(in)    :: ndimEV,numk,ncount
 integer(dp), intent(in)    :: npointsBack(numk,numk)
 real(dp)   , intent(in)    :: fmu,vq1(2), vq12(2), bandsTB(ndimEV,ncount)

 integer(dp) :: i,ivk,ivk1,ivk2
 real(dp) :: vK1x,vK1y,vK2x,vK2y,vM1x,vM1y,vM2x,vM2y
 real(dp) :: a11,a12,a21,a22,b11,b12,b21,b22,vkx,vky,det
 real(dp) :: aGM, aKG, aMK, aT

 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 vK1x=(vq1(1)+vq12(1))/3.0_dp
 vK1y=(vq1(2)+vq12(2))/3.0_dp

 vK2x=2.0_dp*(vq1(1)+vq12(1))/3.0_dp
 vK2y=2.0_dp*(vq1(2)+vq12(2))/3.0_dp

 vM1x=vq1(1)/2.0_dp
 vM1y=vq1(2)/2.0_dp

 vM2x=vq12(1)/2.0_dp
 vM2y=vq12(2)/2.0_dp


 aKG=sqrt(vK1x**2+vK1y**2)
 aGM=sqrt(vM1x**2+vM1y**2)
 aMK=sqrt((vK1x-vM1x)**2+(vK1y-vM1y)**2)
 aT=aKG+aGM+aMK

 a11=vq1(1)/real(numk,dp)
 a12=vq12(1)/real(numk,dp)
 a21=vq1(2)/real(numk,dp)
 a22=vq12(2)/real(numk,dp)

 det=a11*a22-a12*a21
 b11=a22/det
 b12=-a12/det
 b21=-a21/det
 b22=a11/det

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

 do ivk2=0,numk/3
ivk1=ivk2
    write(99,fmt='(F16.8,3X)', advance="no") ivk2*aKG/real(numk/3,dp)/aT
    do i=1,ndimEV-1
    write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(ivk1+1,ivk2+1))-fmu
    enddo
    write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(ivk1+1,ivk2+1))-fmu
enddo
do ivk=numk/6,0,-1
vkx=vM1x+(vK1x-vM1x)*ivk/real(numk/6,dp)
vky=vM2y+(vK1y-vM1y)*ivk/real(numk/6,dp)
ivk2=nint(b11*vkx+b12*vky)
ivk1=nint(b21*vkx+b22*vky)
    write(99,fmt='(F16.8,3X)', advance="no") (aKG+(numk/6-ivk)*aMK/real(numk/6,dp))/aT
    do i=1,ndimEV-1
    write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(ivk1+1,ivk2+1))-fmu
    enddo
    write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(ivk1+1,ivk2+1))-fmu
enddo
do ivk2=numk/2,0,-1
ivk1=0
write(99,fmt='(F16.8,3X)', advance="no") 1.-ivk2*aGM/real(numk/2,dp)/aT
    do i=1,ndimEV-1
    write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(ivk1+1,ivk2+1))-fmu
    enddo
    write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(ivk1+1,ivk2+1))-fmu
enddo
close(99)

end subroutine plotBandsGamma

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine plotBandsAll(bandsTB,fmu,npointsBack,ndimEV,numk,ncount,vq1,vq12)

 integer(dp), intent(in)    :: ndimEV,numk,ncount
 integer(dp), intent(in)    :: npointsBack(numk,numk)
 real(dp)   , intent(in)    :: fmu,vq1(2), vq12(2), bandsTB(ndimEV,ncount)

 integer(dp) :: i,ivk,ivk1,ivk2
 real(dp) :: vK1x,vK1y,vK2x,vK2y,vM1x,vM1y,vM2x,vM2y
 real(dp) :: a11,a12,a21,a22,b11,b12,b21,b22,vkx,vky,det
 real(dp) :: aGM, aKG, aMK, aT

 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 vK1x=(vq1(1)+vq12(1))/3.0_dp
 vK1y=(vq1(2)+vq12(2))/3.0_dp

 vK2x=2.0_dp*(vq1(1)+vq12(1))/3.0_dp
 vK2y=2.0_dp*(vq1(2)+vq12(2))/3.0_dp

 vM1x=vq1(1)/2.0_dp
 vM1y=vq1(2)/2.0_dp

 vM2x=vq12(1)/2.0_dp+vq1(1)
 vM2y=vq12(2)/2.0_dp+vq1(2)


 aKG=sqrt(vK1x**2+vK1y**2)
 aGM=sqrt(vM1x**2+vM1y**2)
 aMK=sqrt((vK1x-vM1x)**2+(vK1y-vM1y)**2)
 aT=aKG+aGM+aMK


 a11=vq1(1)/real(numk,dp)
 a12=vq12(1)/real(numk,dp)
 a21=vq1(2)/real(numk,dp)
 a22=vq12(2)/real(numk,dp)

 det=a11*a22-a12*a21
 b11=a22/det
 b12=-a12/det
 b21=-a21/det
 b22=a11/det

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! The path crosses a reciprocal-cell edge (index numk). Wrap into the
 ! stored 0:numk-1 grid before indexing; the band energies are periodic.
 do ivk2=0,numk/3
ivk1=ivk2
    !write(*,*) '1',ivk2*aKG/real(numk/3,dp)/aT,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
    write(99,fmt='(F16.8,3X)', advance="no") ivk2*aKG/real(numk/3,dp)/aT
    do i=1,ndimEV-1
    write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(modulo(ivk1,numk)+1,modulo(ivk2,numk)+1))-fmu
    enddo
    write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(modulo(ivk1,numk)+1,modulo(ivk2,numk)+1))-fmu
enddo
do ivk=numk/6,0,-1
vkx=vM1x+(vK1x-vM1x)*ivk/real(numk/6,dp)
vky=vM1y+(vK1y-vM1y)*ivk/real(numk/6,dp)
ivk2=nint(b11*vkx+b12*vky)
ivk1=nint(b21*vkx+b22*vky)
    !write(*,*) '2',(aKG+(numk/6-ivk)*aMK/real(numk/6,dp))/aT,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
    write(99,fmt='(F16.8,3X)', advance="no") (aKG+(numk/6-ivk)*aMK/real(numk/6,dp))/aT
    do i=1,ndimEV-1
    write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(modulo(ivk1,numk)+1,modulo(ivk2,numk)+1))-fmu
    enddo
    write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(modulo(ivk1,numk)+1,modulo(ivk2,numk)+1))-fmu
enddo

do ivk1=1,numk/2-1
ivk2=0
!write(*,*) '3',(aKG+aMK+ivk1*aGM/real(numk/2,dp))/aT,((ivk1+numk/2)*vq1(:)+ivk2*vq12(:))/real(numk)
write(99,fmt='(F16.8,3X)', advance="no") (aKG+aMK+ivk1*aGM/real(numk/2,dp))/aT
    do i=1,ndimEV-1
    write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(modulo(ivk1+numk/2,numk)+1,modulo(ivk2,numk)+1))-fmu
    enddo
    write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(modulo(ivk1+numk/2,numk)+1,modulo(ivk2,numk)+1))-fmu
enddo

do ivk2=0,numk/2
ivk1=numk
!write(*,*) '4',1.+ivk2*aGM/real(numk/2,dp)/aT,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
write(99,fmt='(F16.8,3X)', advance="no") 1.+ivk2*aGM/real(numk/2,dp)/aT
    do i=1,ndimEV-1
    write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(modulo(ivk1,numk)+1,modulo(ivk2,numk)+1))-fmu
    enddo
    write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(modulo(ivk1,numk)+1,modulo(ivk2,numk)+1))-fmu
enddo

do ivk=1,numk/6
vkx=vM2x+(vK2x-vM2x)*ivk/real(numk/6,dp)
vky=vM2y+(vK2y-vM2y)*ivk/real(numk/6,dp)
ivk2=nint(b11*vkx+b12*vky)
ivk1=nint(b21*vkx+b22*vky)
!write(*,*) '5',1.+(aGM+ivk*aMK/real(numk/6,dp))/aT,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
    write(99,fmt='(F16.8,3X)', advance="no") 1.+(aGM+ivk*aMK/real(numk/6,dp))/aT
    do i=1,ndimEV-1
    write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(modulo(ivk1,numk)+1,modulo(ivk2,numk)+1))-fmu
    enddo
    write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(modulo(ivk1,numk)+1,modulo(ivk2,numk)+1))-fmu
enddo
 do ivk2=0,numk/3-1
ivk1=ivk2
!write(*,*) '6',1.+(aGM+aMK+ivk2*aKG/real(numk/3,dp))/aT,((ivk1+2*numk/3)*vq1(:)+(ivk2+2*numk/3)*vq12(:))/real(numk)
    write(99,fmt='(F16.8,3X)', advance="no") 1.+(aGM+aMK+ivk2*aKG/real(numk/3,dp))/aT
    do i=1,ndimEV-1
    write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(modulo(ivk1+2*numk/3,numk)+1,modulo(ivk2+2*numk/3,numk)+1))-fmu
    enddo
    write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(modulo(ivk1+2*numk/3,numk)+1,modulo(ivk2+2*numk/3,numk)+1))-fmu
enddo

ivk2=0
ivk1=ivk2
!write(*,*) 2.,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
    write(99,fmt='(F16.8,3X)', advance="no") 2.
    do i=1,ndimEV-1
    write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(modulo(ivk1,numk)+1,modulo(ivk2,numk)+1))-fmu
    enddo
    write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(modulo(ivk1,numk)+1,modulo(ivk2,numk)+1))-fmu

 close(99)

end subroutine plotBandsAll


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine plotBandsDensity(bandsTB,npointsBZ,ndimEV,numk,nband,ncount,nNP)

 integer(dp), intent(in)    :: ndimEV,numk,nband,ncount,nNP
 integer(dp), intent(in)    :: npointsBZ(ncount,8)
 real(dp), intent(in)    :: bandsTB(ndimEV,ncount)

 integer(dp) :: i,icount,ivk1,ivk2

 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 do icount=1,ncount

     ivk1=npointsBZ(icount,1)
     ivk2=npointsBZ(icount,2)

     write(98,fmt='(F16.8,3X)', advance="no") real(ivk1,dp)
     write(98,fmt='(F16.8,3X)', advance="no") real(ivk2,dp)
     do i=nNP-nband,nNP+nband
     write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
     enddo
     write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)

     if(ivk1.EQ.0.AND.ivk2.EQ.0)then
     write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
     write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
     do i=nNP-nband,nNP+nband
     write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
     enddo
     write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)

     write(98,fmt='(F16.8,3X)', advance="no") real(ivk1,dp)
     write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
     do i=nNP-nband,nNP+nband
     write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
     enddo
     write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)

     write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
     write(98,fmt='(F16.8,3X)', advance="no") real(ivk2,dp)
     do i=nNP-nband,nNP+nband
     write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
     enddo
     write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)

     elseif(ivk1.EQ.0)then
     write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
     write(98,fmt='(F16.8,3X)', advance="no") real(ivk2,dp)
     do i=nNP-nband,nNP+nband
     write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
     enddo
     write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)

     elseif(ivk2.EQ.0.)then
     write(98,fmt='(F16.8,3X)', advance="no") real(ivk1,dp)
     write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
     do i=nNP-nband,nNP+nband
     write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
     enddo
     write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)
     endif


enddo

end subroutine plotBandsDensity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine plotBandsDensity3(bandsTB,npointsBZ,ndimEV,numk,nband,ncount,nNP)

 integer(dp), intent(in)    :: ndimEV,numk,nband,ncount,nNP
 integer(dp), intent(in)    :: npointsBZ(ncount,8)
 real(dp), intent(in)    :: bandsTB(ndimEV,ncount)

 integer(dp) :: i,icount,ivk1,ivk2

 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 do icount=1,ncount

     ivk1=npointsBZ(icount,1)
     ivk2=npointsBZ(icount,2)

     write(98,fmt='(F16.8,3X)', advance="no") real(ivk1,dp)
     write(98,fmt='(F16.8,3X)', advance="no") real(ivk2,dp)
     do i=nNP-nband,nNP+nband
     write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
     enddo
     write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)

     if(ivk1.EQ.0.AND.ivk2.EQ.0)then
     write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
     write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
     do i=nNP-nband,nNP+nband
     write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
     enddo
     write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)

     write(98,fmt='(F16.8,3X)', advance="no") real(ivk1,dp)
     write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
     do i=nNP-nband,nNP+nband
     write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
     enddo
     write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)

     write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
     write(98,fmt='(F16.8,3X)', advance="no") real(ivk2,dp)
     do i=nNP-nband,nNP+nband
     write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
     enddo
     write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)

     elseif(ivk1.EQ.0)then
     write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
     write(98,fmt='(F16.8,3X)', advance="no") real(ivk2,dp)
     do i=nNP-nband,nNP+nband
     write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
     enddo
     write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)

     elseif(ivk2.EQ.0.)then
     write(98,fmt='(F16.8,3X)', advance="no") real(ivk1,dp)
     write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
     do i=nNP-nband,nNP+nband
     write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
     enddo
     write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)
     endif



if(npointsBZ(icount,4).EQ.3)then

     ivk1=npointsBZ(icount,5)
     ivk2=npointsBZ(icount,6)


         write(98,fmt='(F16.8,3X)', advance="no") real(ivk1,dp)
         write(98,fmt='(F16.8,3X)', advance="no") real(ivk2,dp)
         do i=nNP-nband,nNP+nband
         write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
         enddo
         write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)

         if(ivk1.EQ.0)then
         write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
         write(98,fmt='(F16.8,3X)', advance="no") real(ivk2,dp)
         do i=nNP-nband,nNP+nband
         write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
         enddo
         write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)

         elseif(ivk2.EQ.0.)then
         write(98,fmt='(F16.8,3X)', advance="no") real(ivk1,dp)
         write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
         do i=nNP-nband,nNP+nband
         write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
         enddo
         write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)
         endif





     ivk1=npointsBZ(icount,7)
     ivk2=npointsBZ(icount,8)


          write(98,fmt='(F16.8,3X)', advance="no") real(ivk1,dp)
          write(98,fmt='(F16.8,3X)', advance="no") real(ivk2,dp)
          do i=nNP-nband,nNP+nband
          write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
          enddo
          write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)

          if(ivk1.EQ.0)then
          write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
          write(98,fmt='(F16.8,3X)', advance="no") real(ivk2,dp)
          do i=nNP-nband,nNP+nband
          write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
          enddo
          write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)

          elseif(ivk2.EQ.0.)then
          write(98,fmt='(F16.8,3X)', advance="no") real(ivk1,dp)
          write(98,fmt='(F16.8,3X)', advance="no") real(numk,dp)
          do i=nNP-nband,nNP+nband
          write(98,fmt='(F16.8,3X)', advance="no") bandsTB(i,icount)
          enddo
          write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband+1,icount)
          endif


endif

enddo

end subroutine plotBandsDensity3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine getDrude(bandsTB,fmu,npointsBack,ndimEV,numk,ncount,vq1,vq12,nxy)

 integer(dp), intent(in)    :: ndimEV,numk,ncount,nxy
 integer(dp), intent(in)    :: npointsBack(numk,numk)
 real(dp)   , intent(in)    :: fmu,vq1(2),vq12(2),bandsTB(ndimEV,ncount)

 real(dp) :: OmegaMin,OmegaMax
 integer(dp), parameter :: nspacing=1000
 integer(dp) :: nomega
 integer(dp) :: nshift

 real(dp), allocatable    :: drude(:),temp1(:),PivotEnergy(:,:),PivotMomentum(:,:)


 integer(dp) :: n,i,ivk,ivk1,ivk2,nbounds(4)
 real(dp) :: Ac,sigma0,dk(2)
 real(dp) :: bandsFull(ndimEV,numk+1,numk+1)
 real(dp) :: vEnergy(ndimEV,numk+1,2),vMomentum(2,numk+1,2)


  OmegaMin = 1.200_dp
  OmegaMax = -1.200_dp

!  do i=1,ncount
!  if(maxval(bandsTB(ndimEV/2-1:ndimEV/2+2,i))-fmu.GT.OmegaMax)then
!  OmegaMax=maxval(bandsTB(ndimEV/2-1:ndimEV/2+2,i))-fmu
!  endif
!  if(minval(bandsTB(ndimEV/2-1:ndimEV/2+2,i))-fmu.LT.OmegaMin)then
!  OmegaMin=minval(bandsTB(ndimEV/2-1:ndimEV/2+2,i))-fmu
!  endif
!  enddo

  do i=1,ncount
  if(maxval(bandsTB(:,i))-fmu.GT.OmegaMax)then
  OmegaMax=maxval(bandsTB(:,i))-fmu
  endif
  if(minval(bandsTB(:,i))-fmu.LT.OmegaMin)then
  OmegaMin=minval(bandsTB(:,i))-fmu
  endif
  enddo

  !write(*,*) 'maxmin',OmegaMax,OmegaMin

  nomega=int((OmegaMax-OmegaMin)*nspacing)
  nshift=int(-OmegaMin*nspacing)


 allocate(drude(nomega),temp1(nomega))

 allocate(pivotEnergy(ndimEV,4),pivotMomentum(2,4))


 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!! Periodic continuation
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 do ivk1=0,numk-1
 do ivk2=0,numk-1
 bandsFull(:,ivk1+1,ivk2+1)=bandsTB(:,npointsBack(ivk1+1,ivk2+1))-fmu
 enddo
 enddo

 do ivk=0,numk-1
 bandsFull(:,ivk+1,numk+1)=bandsFull(:,ivk+1,1)
 bandsFull(:,numk+1,ivk+1)=bandsFull(:,1,ivk+1)
 enddo
 bandsFull(:,numk+1,numk+1)=bandsFull(:,1,1)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

drude=0._dp

nbounds(1)=1 ! min for i_v
nbounds(2)=ndimEV-ncb ! max for i_v
nbounds(3)=ndimEV-ncb+1 ! min for i_c
nbounds(4)=ndimEV ! max for i_c


dk(1)=norm2(vq1(:)+vq12(:))/real(numk,dp)
dk(2)=norm2(vq1(:)-vq12(:))/real(numk,dp)


do ivk1=0,numk-1
 do ivk2 = 0,numk
   vEnergy(:,ivk2+1,1)=bandsFull(:,ivk1+1,ivk2+1)
   vEnergy(:,ivk2+1,2)=bandsFull(:,ivk1+2,ivk2+1)
   vMomentum(:,ivk2+1,1)=(ivk1*vq1(:)+ivk2*vq12(:))/real(numk,dp)
   vMomentum(:,ivk2+1,2)=((ivk1+1)*vq1(:)+ivk2*vq12(:))/real(numk,dp)
 end do

 temp1=0._dp


!$omp parallel do &
!$omp private(ivk2,pivotEnergy,pivotMomentum,temp1) &
!$omp shared(vEnergy,vMomentum,nshift,nomega,nbounds,nxy,dk) &
!$omp reduction(+:drude)
do ivk2=0,numk-1

pivotEnergy=reshape([vEnergy(:,ivk2+1,1),vEnergy(:,ivk2+1,2),vEnergy(:,ivk2+2,1),vEnergy(:,ivk2+2,2)],[ndimEV,4_dp])

pivotMomentum=reshape([vMomentum(:,ivk2+1,1),vMomentum(:,ivk2+1,2),vMomentum(:,ivk2+2,1),vMomentum(:,ivk2+2,2)],[2_dp,4_dp])

call drude_S(temp1,pivotMomentum,pivotEnergy,nomega,nspacing,nshift,nbounds,nxy,dk)

!call dosVB_S(ndim,temp2(:),pivotMomentum,pivotEnergy, &
!nomega,nspacing,nshift,fmu,nbounds)

!call dosJoint_S(ndim,temp3(:),pivotMomentum,pivotEnergy, &
!nomega,nspacing,nshift,nbounds)

drude(:)=drude(:)+temp1(:)
!dosVB(:)=dosVB(:)+temp2(:)
!dosJoint(:)=dosJoint(:)+temp3(:)


enddo ! ivk2
!$omp end parallel do

enddo ! ivk1



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

Ac=4._dp*pi*pi ! normalization: int DOS*sqrt(3)/2*ndim/4= number of bands =iu-il+1



        !do n=nomega,nshift,-1

           !write(99,*) -real(n-nshift)/real(nspacing),dosVB(n)/Ac

        !end do

        do n=1,nomega

           write(99,*) real(n-nshift)/real(nspacing),drude(n)/Ac

        end do

        close(99)


end subroutine getDrude

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine getDrudeV(bandsTB,vxTB,vyTB,fmu,npointsBack,ndimEV,numk,ncount,vq1,vq12,coord,tn)

 integer(dp), intent(in)    :: ndimEV,numk,ncount
 integer(dp), intent(in)    :: npointsBack(numk,numk)
 real(dp)   , intent(in)    :: fmu,vq1(2),vq12(2),bandsTB(ndimEV,ncount)
 real(dp) , intent(in) :: tn(2,3), coord(ndim,3)
 real(dp)   , intent(in)    :: vxTB(ndimEV,ncount),vyTB(ndimEV,ncount)

 real(dp) :: OmegaMin,OmegaMax
 integer(dp), parameter :: nspacing=1000
 integer(dp) :: nomega
 integer(dp) :: nshift

 real(dp), allocatable    :: drudeX(:),drudeY(:),temp1(:),temp2(:),PivotEnergy(:,:),PivotMomentum(:,:),PivotVx(:,:),PivotVy(:,:)


 integer(dp) :: n,i,ivk,ivk1,ivk2,nbounds(4)
 real(dp) :: Ac,sigma0,dk(2)
 real(dp) :: bandsFull(ndimEV,numk+1,numk+1)
 real(dp) :: vEnergy(ndimEV,numk+1,2),vMomentum(2,numk+1,2),vVx(ndimEV,numk+1,2),vVy(ndimEV,numk+1,2)
 real(dp) :: rout(ndimEV,2)
 real(dp) :: vxTBFull(ndimEV,numk+1,numk+1),vyTBFull(ndimEV,numk+1,numk+1)


  OmegaMin = 1.20_dp
  OmegaMax = -1.20_dp

!  do i=1,ncount
!  if(maxval(bandsTB(ndimEV/2-1:ndimEV/2+2,i))-fmu.GT.OmegaMax)then
!  OmegaMax=maxval(bandsTB(ndimEV/2-1:ndimEV/2+2,i))-fmu
!  endif
!  if(minval(bandsTB(ndimEV/2-1:ndimEV/2+2,i))-fmu.LT.OmegaMin)then
!  OmegaMin=minval(bandsTB(ndimEV/2-1:ndimEV/2+2,i))-fmu
!  endif
!  enddo

  do i=1,ncount
  if(maxval(bandsTB(:,i))-fmu.GT.OmegaMax)then
  OmegaMax=maxval(bandsTB(:,i))-fmu
  endif
  if(minval(bandsTB(:,i))-fmu.LT.OmegaMin)then
  OmegaMin=minval(bandsTB(:,i))-fmu
  endif
  enddo

  write(*,*) 'maxmin',OmegaMax,OmegaMin

  nomega=int((OmegaMax-OmegaMin)*nspacing)
  nshift=int(-OmegaMin*nspacing)


 allocate(drudeX(nomega),drudeY(nomega),temp1(nomega),temp2(nomega))

 allocate(pivotEnergy(ndimEV,4),pivotMomentum(2,4),pivotVx(ndimEV,4),pivotVy(ndimEV,4))


 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!! Periodic continuation
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 do ivk1=0,numk-1
 do ivk2=0,numk-1
 bandsFull(:,ivk1+1,ivk2+1)=bandsTB(:,npointsBack(ivk1+1,ivk2+1))-fmu
 vxTBFull(:,ivk1+1,ivk2+1)=vxTB(:,npointsBack(ivk1+1,ivk2+1))
 vyTBFull(:,ivk1+1,ivk2+1)=vyTB(:,npointsBack(ivk1+1,ivk2+1))
 enddo
 enddo

 do ivk=0,numk-1
 bandsFull(:,ivk+1,numk+1)=bandsFull(:,ivk+1,1)
 bandsFull(:,numk+1,ivk+1)=bandsFull(:,1,ivk+1)
 vxTBFull(:,ivk+1,numk+1)=vxTBFull(:,ivk+1,1)
 vxTBFull(:,numk+1,ivk+1)=vxTBFull(:,1,ivk+1)
 vyTBFull(:,ivk+1,numk+1)=vyTBFull(:,ivk+1,1)
 vyTBFull(:,numk+1,ivk+1)=vyTBFull(:,1,ivk+1)
 enddo
 bandsFull(:,numk+1,numk+1)=bandsFull(:,1,1)
 vxTBFull(:,numk+1,numk+1)=vxTBFull(:,1,1)
 vyTBFull(:,numk+1,numk+1)=vyTBFull(:,1,1)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

drudeX=0._dp
drudeY=0._dp

nbounds(1)=1 ! min for i_v
nbounds(2)=ndimEV-ncb ! max for i_v
nbounds(3)=ndimEV-ncb+1 ! min for i_c
nbounds(4)=ndimEV ! max for i_c

ivk1=0
 do ivk2 = 0,numk
 write(*,*) ivk2
   vEnergy(:,ivk2+1,1)=bandsFull(:,ivk1+1,ivk2+1)
   vMomentum(:,ivk2+1,1)=(ivk1*vq1(:)+ivk2*vq12(:))/real(numk,dp)
   !call overlap(rout,zevFull(:,:,ivk1+1,ivk2+1),coord,vMomentum(:,ivk2+1,1),tn,ndimEV)
   vVx(:,ivk2+1,1)=vxTBFull(:,ivk1+1,ivk2+1)
   vVy(:,ivk2+1,1)=vyTBFull(:,ivk1+1,ivk2+1)
 end do

do ivk1=1,numk
 write(*,*) ivk1
do ivk2 = 0,numk
   vEnergy(:,ivk2+1,2)=bandsFull(:,ivk1+1,ivk2+1)
   vMomentum(:,ivk2+1,2)=(ivk1*vq1(:)+ivk2*vq12(:))/real(numk,dp)
   !call overlap(rout,zevFull(:,:,ivk1+1,ivk2+1),coord,vMomentum(:,ivk2+1,1),tn,ndimEV)
   vVx(:,ivk2+1,2)=vxTBFull(:,ivk1+1,ivk2+1)
   vVy(:,ivk2+1,2)=vyTBFull(:,ivk1+1,ivk2+1)
 end do

 temp1=0._dp
 temp2=0._dp

!$omp parallel do &
!$omp private(ivk2,pivotEnergy,pivotMomentum,temp1,temp2) &
!$omp shared(vEnergy,vMomentum,vVx,vVy,nshift,nomega,nbounds) &
!$omp reduction(+:drudeX) &
!$omp reduction(+:drudeY)
do ivk2=0,numk-1

pivotEnergy=reshape([vEnergy(:,ivk2+1,1),vEnergy(:,ivk2+1,2),vEnergy(:,ivk2+2,1),vEnergy(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVx=reshape([vVx(:,ivk2+1,1),vVx(:,ivk2+1,2),vVx(:,ivk2+2,1),vVx(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVy=reshape([vVy(:,ivk2+1,1),vVy(:,ivk2+1,2),vVy(:,ivk2+2,1),vVy(:,ivk2+2,2)],[ndimEV,4_dp])

pivotMomentum=reshape([vMomentum(:,ivk2+1,1),vMomentum(:,ivk2+1,2),vMomentum(:,ivk2+2,1),vMomentum(:,ivk2+2,2)],[2_dp,4_dp])

call drude_V(temp1,temp2,pivotMomentum,pivotEnergy,pivotVx,pivotVy,nomega,nspacing,nshift,nbounds)

!call dosVB_S(ndim,temp2(:),pivotMomentum,pivotEnergy, &
!nomega,nspacing,nshift,fmu,nbounds)

!call dosJoint_S(ndim,temp3(:),pivotMomentum,pivotEnergy, &
!nomega,nspacing,nshift,nbounds)

drudeX(:)=drudeX(:)+temp1(:)
drudeY(:)=drudeY(:)+temp2(:)
!dosVB(:)=dosVB(:)+temp2(:)
!dosJoint(:)=dosJoint(:)+temp3(:)
enddo ! ivk2
!$omp end parallel do
vEnergy(:,:,1)=vEnergy(:,:,2)
vMomentum(:,:,1)=vMomentum(:,:,2)
vVx(:,:,1)=vVx(:,:,2)
vVy(:,:,1)=vVy(:,:,2)
write(*,*) ivk1
enddo ! ivk1



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

Ac=4._dp*pi*pi ! normalization: int DOS*sqrt(3)/2*ndim/4= number of bands =iu-il+1


write(*,*) 'TestWriteDrude'
        !do n=nomega,nshift,-1

           !write(99,*) -real(n-nshift)/real(nspacing),dosVB(n)/Ac

        !end do

        do n=1,nomega

           write(98,*) real(n-nshift)/real(nspacing),drudeX(n)/Ac
           write(99,*) real(n-nshift)/real(nspacing),drudeY(n)/Ac

        end do

        close(98)
        close(99)


end subroutine getDrudeV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine getDrudeFullvc(bandsTB,drudeTotVB,drudeMagVB,drudeChiVB,drudeTotCB,drudeMagCB,drudeChiCB,vx1TB,vy1TB,vx2TB,vy2TB,vxITB,vyITB, &
fmu,npointsBack,ndimEV,numk,ncount,vq1,vq12,nomega,nspacing,nshiftVB,nshiftCB,nbounds)

 integer(dp), intent(in)    :: ndimEV,numk,ncount
 integer(dp), intent(in)    :: nomega,nspacing,nshiftVB,nshiftCB,nbounds(:)
 integer(dp), intent(in)    :: npointsBack(numk,numk)
 real(dp)   , intent(in)    :: fmu,vq1(2),vq12(2),bandsTB(ndimEV,ncount)
 real(dp)   , intent(in)    :: vx1TB(:,:),vy1TB(:,:),vx2TB(:,:),vy2TB(:,:),vxITB(:,:),vyITB(:,:)
 real(dp)   , intent(out)    :: drudeTotVB(:),drudeMagVB(:),drudeChiVB(:)
 real(dp)   , intent(out)    :: drudeTotCB(:),drudeMagCB(:),drudeChiCB(:)


 real(dp), allocatable    :: temp1(:),temp2(:),temp3(:)
 real(dp), allocatable    :: PivotEnergy(:,:),PivotMomentum(:,:),PivotVx1(:,:),PivotVy1(:,:),PivotVx2(:,:),PivotVy2(:,:),PivotVxI(:,:),PivotVyI(:,:)


 integer(dp) :: n,i,ivk,ivk1,ivk2
 real(dp) :: Ac,sigma0,dk(2)
 real(dp) :: bandsFull(ndimEV,numk+1,numk+1)
 real(dp) :: vEnergy(ndimEV,numk+1,2),vMomentum(2,numk+1,2)
 real(dp) :: vVx1(ndimEV,numk+1,2),vVy1(ndimEV,numk+1,2),vVx2(ndimEV,numk+1,2),vVy2(ndimEV,numk+1,2),vVxI(ndimEV,numk+1,2),vVyI(ndimEV,numk+1,2)
 real(dp) :: rout(ndimEV,2)
 real(dp) :: vx1TBFull(ndimEV,numk+1,numk+1),vy1TBFull(ndimEV,numk+1,numk+1)
 real(dp) :: vx2TBFull(ndimEV,numk+1,numk+1),vy2TBFull(ndimEV,numk+1,numk+1)
 real(dp) :: vxITBFull(ndimEV,numk+1,numk+1),vyITBFull(ndimEV,numk+1,numk+1)






 allocate(temp1(nomega),temp2(nomega),temp3(nomega))

 allocate(pivotEnergy(ndimEV,4),pivotMomentum(2,4))
 allocate(pivotVx1(ndimEV,4),pivotVy1(ndimEV,4),pivotVx2(ndimEV,4),pivotVy2(ndimEV,4),pivotVxI(ndimEV,4),pivotVyI(ndimEV,4))


 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!! Periodic continuation
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 do ivk1=0,numk-1
 do ivk2=0,numk-1
 bandsFull(:,ivk1+1,ivk2+1)=bandsTB(:,npointsBack(ivk1+1,ivk2+1))
 vx1TBFull(:,ivk1+1,ivk2+1)=vx1TB(:,npointsBack(ivk1+1,ivk2+1))
 vy1TBFull(:,ivk1+1,ivk2+1)=vy1TB(:,npointsBack(ivk1+1,ivk2+1))
 vx2TBFull(:,ivk1+1,ivk2+1)=vx2TB(:,npointsBack(ivk1+1,ivk2+1))
 vy2TBFull(:,ivk1+1,ivk2+1)=vy2TB(:,npointsBack(ivk1+1,ivk2+1))
 vxITBFull(:,ivk1+1,ivk2+1)=vxITB(:,npointsBack(ivk1+1,ivk2+1))
 vyITBFull(:,ivk1+1,ivk2+1)=vyITB(:,npointsBack(ivk1+1,ivk2+1))
 enddo
 enddo

 do ivk=0,numk-1
 bandsFull(:,ivk+1,numk+1)=bandsFull(:,ivk+1,1)
 bandsFull(:,numk+1,ivk+1)=bandsFull(:,1,ivk+1)
 vx1TBFull(:,ivk+1,numk+1)=vx1TBFull(:,ivk+1,1)
 vx1TBFull(:,numk+1,ivk+1)=vx1TBFull(:,1,ivk+1)
 vy1TBFull(:,ivk+1,numk+1)=vy1TBFull(:,ivk+1,1)
 vy1TBFull(:,numk+1,ivk+1)=vy1TBFull(:,1,ivk+1)
 vx2TBFull(:,ivk+1,numk+1)=vx2TBFull(:,ivk+1,1)
 vx2TBFull(:,numk+1,ivk+1)=vx2TBFull(:,1,ivk+1)
 vy2TBFull(:,ivk+1,numk+1)=vy2TBFull(:,ivk+1,1)
 vy2TBFull(:,numk+1,ivk+1)=vy2TBFull(:,1,ivk+1)
 vxITBFull(:,ivk+1,numk+1)=vxITBFull(:,ivk+1,1)
 vxITBFull(:,numk+1,ivk+1)=vxITBFull(:,1,ivk+1)
 vyITBFull(:,ivk+1,numk+1)=vyITBFull(:,ivk+1,1)
 vyITBFull(:,numk+1,ivk+1)=vyITBFull(:,1,ivk+1)
 enddo
 bandsFull(:,numk+1,numk+1)=bandsFull(:,1,1)
 vx1TBFull(:,numk+1,numk+1)=vx1TBFull(:,1,1)
 vy1TBFull(:,numk+1,numk+1)=vy1TBFull(:,1,1)
 vx2TBFull(:,numk+1,numk+1)=vx2TBFull(:,1,1)
 vy2TBFull(:,numk+1,numk+1)=vy2TBFull(:,1,1)
 vxITBFull(:,numk+1,numk+1)=vxITBFull(:,1,1)
 vyITBFull(:,numk+1,numk+1)=vyITBFull(:,1,1)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

drudeTotVB=0._dp
drudeMagVB=0._dp
drudeChiVB=0._dp

drudeTotCB=0._dp
drudeMagCB=0._dp
drudeChiCB=0._dp


ivk1=0
 do ivk2 = 0,numk
 write(*,*) ivk2
   vEnergy(:,ivk2+1,1)=bandsFull(:,ivk1+1,ivk2+1)
   vMomentum(:,ivk2+1,1)=(ivk1*vq1(:)+ivk2*vq12(:))/real(numk,dp)
   vVx1(:,ivk2+1,1)=vx1TBFull(:,ivk1+1,ivk2+1)
   vVy1(:,ivk2+1,1)=vy1TBFull(:,ivk1+1,ivk2+1)
   vVx2(:,ivk2+1,1)=vx2TBFull(:,ivk1+1,ivk2+1)
   vVy2(:,ivk2+1,1)=vy2TBFull(:,ivk1+1,ivk2+1)
   vVxI(:,ivk2+1,1)=vxITBFull(:,ivk1+1,ivk2+1)
   vVyI(:,ivk2+1,1)=vyITBFull(:,ivk1+1,ivk2+1)
 end do

do ivk1=1,numk
 write(*,*) ivk1
do ivk2 = 0,numk
   vEnergy(:,ivk2+1,2)=bandsFull(:,ivk1+1,ivk2+1)
   vMomentum(:,ivk2+1,2)=(ivk1*vq1(:)+ivk2*vq12(:))/real(numk,dp)
   vVx1(:,ivk2+1,2)=vx1TBFull(:,ivk1+1,ivk2+1)
   vVy1(:,ivk2+1,2)=vy1TBFull(:,ivk1+1,ivk2+1)
   vVx2(:,ivk2+1,2)=vx2TBFull(:,ivk1+1,ivk2+1)
   vVy2(:,ivk2+1,2)=vy2TBFull(:,ivk1+1,ivk2+1)
   vVxI(:,ivk2+1,2)=vxITBFull(:,ivk1+1,ivk2+1)
   vVyI(:,ivk2+1,2)=vyITBFull(:,ivk1+1,ivk2+1)
 end do


!$omp parallel do &
!$omp private(ivk2,pivotEnergy,pivotMomentum,temp1,temp2) &
!$omp private(pivotVx1,pivotVy1,pivotVx2,pivotVy2,pivotVxI,pivotVyI) &
!$omp shared(vEnergy,vMomentum,vVx1,vVy1,vVx2,vVy2,vVxI,vVyI,nshiftVB,nshiftCB,nomega,nbounds) &
!$omp reduction(+:drudeTotVB) &
!$omp reduction(+:drudeMagVB) &
!$omp reduction(+:drudeChiVB) &
!$omp reduction(+:drudeTotCB) &
!$omp reduction(+:drudeMagCB) &
!$omp reduction(+:drudeChiCB)
do ivk2=0,numk-1

pivotEnergy=reshape([vEnergy(:,ivk2+1,1),vEnergy(:,ivk2+1,2),vEnergy(:,ivk2+2,1),vEnergy(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVx1=reshape([vVx1(:,ivk2+1,1),vVx1(:,ivk2+1,2),vVx1(:,ivk2+2,1),vVx1(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVy1=reshape([vVy1(:,ivk2+1,1),vVy1(:,ivk2+1,2),vVy1(:,ivk2+2,1),vVy1(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVx2=reshape([vVx2(:,ivk2+1,1),vVx2(:,ivk2+1,2),vVx2(:,ivk2+2,1),vVx2(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVy2=reshape([vVy2(:,ivk2+1,1),vVy2(:,ivk2+1,2),vVy2(:,ivk2+2,1),vVy2(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVxI=reshape([vVxI(:,ivk2+1,1),vVxI(:,ivk2+1,2),vVxI(:,ivk2+2,1),vVxI(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVyI=reshape([vVyI(:,ivk2+1,1),vVyI(:,ivk2+1,2),vVyI(:,ivk2+2,1),vVyI(:,ivk2+2,2)],[ndimEV,4_dp])

pivotMomentum=reshape([vMomentum(:,ivk2+1,1),vMomentum(:,ivk2+1,2),vMomentum(:,ivk2+2,1),vMomentum(:,ivk2+2,2)],[2_dp,4_dp])

call drudeVB(temp1,temp2,temp3,pivotMomentum,pivotEnergy,pivotVx1,pivotVy1,pivotVx2,pivotVy2,pivotVxI,pivotVyI,nomega,nspacing,nshiftVB,fmu,nbounds)

drudeTotVB(:)=drudeTotVB(:)+temp1(:)
drudeMagVB(:)=drudeMagVB(:)+temp2(:)
drudeChiVB(:)=drudeChiVB(:)+temp3(:)

call drudeCB(temp1,temp2,temp3,pivotMomentum,pivotEnergy,pivotVx1,pivotVy1,pivotVx2,pivotVy2,pivotVxI,pivotVyI,nomega,nspacing,nshiftCB,fmu,nbounds)

drudeTotCB(:)=drudeTotCB(:)+temp1(:)
drudeMagCB(:)=drudeMagCB(:)+temp2(:)
drudeChiCB(:)=drudeChiCB(:)+temp3(:)

enddo ! ivk2
!$omp end parallel do
vEnergy(:,:,1)=vEnergy(:,:,2)
vMomentum(:,:,1)=vMomentum(:,:,2)
vVx1(:,:,1)=vVx1(:,:,2)
vVy1(:,:,1)=vVy1(:,:,2)
vVx2(:,:,1)=vVx2(:,:,2)
vVy2(:,:,1)=vVy2(:,:,2)
vVxI(:,:,1)=vVxI(:,:,2)
vVyI(:,:,1)=vVyI(:,:,2)
write(*,*) ivk1
enddo ! ivk1



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

Ac=4._dp*pi*pi ! normalization: int DOS*sqrt(3)/2*ndim/4= number of bands =iu-il+1

        do n=nomega,nshiftVB,-1

           write(99,*) -real(n-nshiftVB)/real(nspacing),drudeTotVB(n)/Ac

        end do

        do n=nshiftCB,nomega

           write(99,*) real(n-nshiftCB)/real(nspacing),drudeTotCB(n)/Ac

        end do

        close(99)



end subroutine getDrudeFullvc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine getDrudeFull(bandsTB,drudeTot,drudeMag,drudeChi,vx1TB,vy1TB,vx2TB,vy2TB,vxITB,vyITB, &
fmu,npointsBack,ndimEV,numk,ncount,vq1,vq12,nomega,nspacing,nshift,nbounds)

 integer(dp), intent(in)    :: ndimEV,numk,ncount
 integer(dp), intent(in)    :: nomega,nspacing,nshift,nbounds(:)
 integer(dp), intent(in)    :: npointsBack(numk,numk)
 real(dp)   , intent(in)    :: fmu,vq1(2),vq12(2),bandsTB(ndimEV,ncount)
 real(dp)   , intent(in)    :: vx1TB(:,:),vy1TB(:,:),vx2TB(:,:),vy2TB(:,:),vxITB(:,:),vyITB(:,:)
 real(dp)   , intent(out)    :: drudeTot(:),drudeMag(:),drudeChi(:)


 real(dp), allocatable    :: temp1(:),temp2(:),temp3(:)
 real(dp), allocatable    :: PivotEnergy(:,:),PivotMomentum(:,:),PivotVx1(:,:),PivotVy1(:,:),PivotVx2(:,:),PivotVy2(:,:),PivotVxI(:,:),PivotVyI(:,:)


 integer(dp) :: n,i,ivk,ivk1,ivk2
 real(dp) :: Ac,sigma0,dk(2)
 real(dp) :: bandsFull(ndimEV,numk+1,numk+1)
 real(dp) :: vEnergy(ndimEV,numk+1,2),vMomentum(2,numk+1,2)
 real(dp) :: vVx1(ndimEV,numk+1,2),vVy1(ndimEV,numk+1,2),vVx2(ndimEV,numk+1,2),vVy2(ndimEV,numk+1,2),vVxI(ndimEV,numk+1,2),vVyI(ndimEV,numk+1,2)
 real(dp) :: rout(ndimEV,2)
 real(dp) :: vx1TBFull(ndimEV,numk+1,numk+1),vy1TBFull(ndimEV,numk+1,numk+1)
 real(dp) :: vx2TBFull(ndimEV,numk+1,numk+1),vy2TBFull(ndimEV,numk+1,numk+1)
 real(dp) :: vxITBFull(ndimEV,numk+1,numk+1),vyITBFull(ndimEV,numk+1,numk+1)






 allocate(temp1(nomega),temp2(nomega),temp3(nomega))

 allocate(pivotEnergy(ndimEV,4),pivotMomentum(2,4))
 allocate(pivotVx1(ndimEV,4),pivotVy1(ndimEV,4),pivotVx2(ndimEV,4),pivotVy2(ndimEV,4),pivotVxI(ndimEV,4),pivotVyI(ndimEV,4))


 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!! Periodic continuation
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 do ivk1=0,numk-1
 do ivk2=0,numk-1
 bandsFull(:,ivk1+1,ivk2+1)=bandsTB(:,npointsBack(ivk1+1,ivk2+1))
 vx1TBFull(:,ivk1+1,ivk2+1)=vx1TB(:,npointsBack(ivk1+1,ivk2+1))
 vy1TBFull(:,ivk1+1,ivk2+1)=vy1TB(:,npointsBack(ivk1+1,ivk2+1))
 vx2TBFull(:,ivk1+1,ivk2+1)=vx2TB(:,npointsBack(ivk1+1,ivk2+1))
 vy2TBFull(:,ivk1+1,ivk2+1)=vy2TB(:,npointsBack(ivk1+1,ivk2+1))
 vxITBFull(:,ivk1+1,ivk2+1)=vxITB(:,npointsBack(ivk1+1,ivk2+1))
 vyITBFull(:,ivk1+1,ivk2+1)=vyITB(:,npointsBack(ivk1+1,ivk2+1))
 enddo
 enddo

 do ivk=0,numk-1
 bandsFull(:,ivk+1,numk+1)=bandsFull(:,ivk+1,1)
 bandsFull(:,numk+1,ivk+1)=bandsFull(:,1,ivk+1)
 vx1TBFull(:,ivk+1,numk+1)=vx1TBFull(:,ivk+1,1)
 vx1TBFull(:,numk+1,ivk+1)=vx1TBFull(:,1,ivk+1)
 vy1TBFull(:,ivk+1,numk+1)=vy1TBFull(:,ivk+1,1)
 vy1TBFull(:,numk+1,ivk+1)=vy1TBFull(:,1,ivk+1)
 vx2TBFull(:,ivk+1,numk+1)=vx2TBFull(:,ivk+1,1)
 vx2TBFull(:,numk+1,ivk+1)=vx2TBFull(:,1,ivk+1)
 vy2TBFull(:,ivk+1,numk+1)=vy2TBFull(:,ivk+1,1)
 vy2TBFull(:,numk+1,ivk+1)=vy2TBFull(:,1,ivk+1)
 vxITBFull(:,ivk+1,numk+1)=vxITBFull(:,ivk+1,1)
 vxITBFull(:,numk+1,ivk+1)=vxITBFull(:,1,ivk+1)
 vyITBFull(:,ivk+1,numk+1)=vyITBFull(:,ivk+1,1)
 vyITBFull(:,numk+1,ivk+1)=vyITBFull(:,1,ivk+1)
 enddo
 bandsFull(:,numk+1,numk+1)=bandsFull(:,1,1)
 vx1TBFull(:,numk+1,numk+1)=vx1TBFull(:,1,1)
 vy1TBFull(:,numk+1,numk+1)=vy1TBFull(:,1,1)
 vx2TBFull(:,numk+1,numk+1)=vx2TBFull(:,1,1)
 vy2TBFull(:,numk+1,numk+1)=vy2TBFull(:,1,1)
 vxITBFull(:,numk+1,numk+1)=vxITBFull(:,1,1)
 vyITBFull(:,numk+1,numk+1)=vyITBFull(:,1,1)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

drudeTot=0._dp
drudeMag=0._dp
drudeChi=0._dp



ivk1=0
 do ivk2 = 0,numk
   vEnergy(:,ivk2+1,1)=bandsFull(:,ivk1+1,ivk2+1)
   vMomentum(:,ivk2+1,1)=(ivk1*vq1(:)+ivk2*vq12(:))/real(numk,dp)
   vVx1(:,ivk2+1,1)=vx1TBFull(:,ivk1+1,ivk2+1)
   vVy1(:,ivk2+1,1)=vy1TBFull(:,ivk1+1,ivk2+1)
   vVx2(:,ivk2+1,1)=vx2TBFull(:,ivk1+1,ivk2+1)
   vVy2(:,ivk2+1,1)=vy2TBFull(:,ivk1+1,ivk2+1)
   vVxI(:,ivk2+1,1)=vxITBFull(:,ivk1+1,ivk2+1)
   vVyI(:,ivk2+1,1)=vyITBFull(:,ivk1+1,ivk2+1)
 end do

do ivk1=1,numk
do ivk2 = 0,numk
   vEnergy(:,ivk2+1,2)=bandsFull(:,ivk1+1,ivk2+1)
   vMomentum(:,ivk2+1,2)=(ivk1*vq1(:)+ivk2*vq12(:))/real(numk,dp)
   vVx1(:,ivk2+1,2)=vx1TBFull(:,ivk1+1,ivk2+1)
   vVy1(:,ivk2+1,2)=vy1TBFull(:,ivk1+1,ivk2+1)
   vVx2(:,ivk2+1,2)=vx2TBFull(:,ivk1+1,ivk2+1)
   vVy2(:,ivk2+1,2)=vy2TBFull(:,ivk1+1,ivk2+1)
   vVxI(:,ivk2+1,2)=vxITBFull(:,ivk1+1,ivk2+1)
   vVyI(:,ivk2+1,2)=vyITBFull(:,ivk1+1,ivk2+1)
 end do


!$omp parallel do &
!$omp private(ivk2,pivotEnergy,pivotMomentum,temp1,temp2) &
!$omp private(pivotVx1,pivotVy1,pivotVx2,pivotVy2,pivotVxI,pivotVyI) &
!$omp shared(vEnergy,vMomentum,vVx1,vVy1,vVx2,vVy2,vVxI,vVyI,nshift,nomega,nbounds) &
!$omp reduction(+:drudeTot) &
!$omp reduction(+:drudeMag) &
!$omp reduction(+:drudeChi)
do ivk2=0,numk-1

pivotEnergy=reshape([vEnergy(:,ivk2+1,1),vEnergy(:,ivk2+1,2),vEnergy(:,ivk2+2,1),vEnergy(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVx1=reshape([vVx1(:,ivk2+1,1),vVx1(:,ivk2+1,2),vVx1(:,ivk2+2,1),vVx1(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVy1=reshape([vVy1(:,ivk2+1,1),vVy1(:,ivk2+1,2),vVy1(:,ivk2+2,1),vVy1(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVx2=reshape([vVx2(:,ivk2+1,1),vVx2(:,ivk2+1,2),vVx2(:,ivk2+2,1),vVx2(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVy2=reshape([vVy2(:,ivk2+1,1),vVy2(:,ivk2+1,2),vVy2(:,ivk2+2,1),vVy2(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVxI=reshape([vVxI(:,ivk2+1,1),vVxI(:,ivk2+1,2),vVxI(:,ivk2+2,1),vVxI(:,ivk2+2,2)],[ndimEV,4_dp])

pivotVyI=reshape([vVyI(:,ivk2+1,1),vVyI(:,ivk2+1,2),vVyI(:,ivk2+2,1),vVyI(:,ivk2+2,2)],[ndimEV,4_dp])

pivotMomentum=reshape([vMomentum(:,ivk2+1,1),vMomentum(:,ivk2+1,2),vMomentum(:,ivk2+2,1),vMomentum(:,ivk2+2,2)],[2_dp,4_dp])

call drudeFull(temp1,temp2,temp3,pivotMomentum,pivotEnergy,pivotVx1,pivotVy1,pivotVx2,pivotVy2,pivotVxI,pivotVyI,nomega,nspacing,nshift,fmu,nbounds)

drudeTot(:)=drudeTot(:)+temp1(:)
drudeMag(:)=drudeMag(:)+temp2(:)
drudeChi(:)=drudeChi(:)+temp3(:)

enddo ! ivk2
!$omp end parallel do
vEnergy(:,:,1)=vEnergy(:,:,2)
vMomentum(:,:,1)=vMomentum(:,:,2)
vVx1(:,:,1)=vVx1(:,:,2)
vVy1(:,:,1)=vVy1(:,:,2)
vVx2(:,:,1)=vVx2(:,:,2)
vVy2(:,:,1)=vVy2(:,:,2)
vVxI(:,:,1)=vVxI(:,:,2)
vVyI(:,:,1)=vVyI(:,:,2)
enddo ! ivk1



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

Ac=4._dp*pi*pi ! normalization: int DOS*sqrt(3)/2*ndim/4= number of bands =iu-il+1

        do n=1,nomega

           write(99,*) real(n-nshift)/real(nspacing),drudeTot(n)/Ac

        end do

        close(99)

end subroutine getDrudeFull



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dosFullvc(bandsTB,dosVB,dosCB,fmu,npointsBack,ndimEV,numk,ncount,vq1,vq12,nomega,nspacing,nshiftVB,nshiftCB,nbounds)

 integer(dp), intent(in)    :: ndimEV,numk,ncount
 integer(dp), intent(in)    :: nomega,nspacing,nshiftVB,nshiftCB,nbounds(:)
 integer(dp), intent(in)    :: npointsBack(numk,numk)
 real (dp)   , intent(in)    :: fmu,vq1(2),vq12(2),bandsTB(ndimEV,ncount)
 real(dp), intent(out) :: dosVB(:),dosCB(:)


 real(dp), allocatable    :: temp1(:),temp2(:),PivotEnergy(:,:),PivotMomentum(:,:)


 integer(dp) :: n,i,ivk,ivk1,ivk2
 real(dp) :: Ac,sigma0
 real(dp) :: bandsFull(ndimEV,numk+1,numk+1)
 real(dp) :: vEnergy(ndimEV,numk+1,2),vMomentum(2,numk+1,2)



!  do i=1,ncount
!  if(maxval(bandsTB(ndimEV/2-1:ndimEV/2+2,i))-fmu.GT.OmegaMax)then
!  OmegaMax=maxval(bandsTB(ndimEV/2-1:ndimEV/2+2,i))-fmu
!  endif
!  if(minval(bandsTB(ndimEV/2-1:ndimEV/2+2,i))-fmu.LT.OmegaMin)then
!  OmegaMin=minval(bandsTB(ndimEV/2-1:ndimEV/2+2,i))-fmu
!  endif
!  enddo


 allocate(temp1(nomega),temp2(nomega))

 allocate(pivotEnergy(ndimEV,4),pivotMomentum(2,4))


 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!! Periodic continuation
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 do ivk1=0,numk-1
 do ivk2=0,numk-1
 bandsFull(:,ivk1+1,ivk2+1)=bandsTB(:,npointsBack(ivk1+1,ivk2+1))
 enddo
 enddo

 do ivk=0,numk-1
 bandsFull(:,ivk+1,numk+1)=bandsFull(:,ivk+1,1)
 bandsFull(:,numk+1,ivk+1)=bandsFull(:,1,ivk+1)
 enddo
 bandsFull(:,numk+1,numk+1)=bandsFull(:,1,1)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

dosCB=0._dp
dosVB=0._dp


do ivk1=0,numk-1
 do ivk2 = 0,numk
   vEnergy(:,ivk2+1,1)=bandsFull(:,ivk1+1,ivk2+1)
   vEnergy(:,ivk2+1,2)=bandsFull(:,ivk1+2,ivk2+1)
   vMomentum(:,ivk2+1,1)=(ivk1*vq1(:)+ivk2*vq12(:))/real(numk,dp)
   vMomentum(:,ivk2+1,2)=((ivk1+1)*vq1(:)+ivk2*vq12(:))/real(numk,dp)
 end do

 temp1=0._dp
 temp2=0._dp

!$omp parallel do &
!$omp private(ivk2,pivotEnergy,pivotMomentum,temp1,temp2) &
!$omp shared(vEnergy,vMomentum,nshiftVB,nshiftCB,nomega,nbounds) &
!$omp reduction(+:dosCB) &
!$omp reduction(+:dosVB)
!!$omp reduction(+:dosJoint)
do ivk2=0,numk-1

pivotEnergy=reshape([vEnergy(:,ivk2+1,1),vEnergy(:,ivk2+1,2),vEnergy(:,ivk2+2,1),vEnergy(:,ivk2+2,2)],[ndimEV,4_dp])

pivotMomentum=reshape([vMomentum(:,ivk2+1,1),vMomentum(:,ivk2+1,2),vMomentum(:,ivk2+2,1),vMomentum(:,ivk2+2,2)],[2_dp,4_dp])

!call dos_S(temp1,pivotMomentum,pivotEnergy,nomega,nspacing,nshift,nbounds)

call dosCB_S(temp1,pivotMomentum,pivotEnergy,nomega,nspacing,nshiftCB,fmu,nbounds)

call dosVB_S(temp2,pivotMomentum,pivotEnergy,nomega,nspacing,nshiftVB,fmu,nbounds)

!call dosJoint_S(ndim,temp3(:),pivotMomentum,pivotEnergy, &
!nomega,nspacing,nshift,nbounds)

dosCB(:)=dosCB(:)+temp1(:)
dosVB(:)=dosVB(:)+temp2(:)
!dosJoint(:)=dosJoint(:)+temp3(:)


enddo ! ivk2
!$omp end parallel do

enddo ! ivk1



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

Ac=4._dp*pi*pi ! normalization: int DOS*sqrt(3)/2*ndim/4= number of bands =iu-il+1


        !do n=nomega,nshift,-1

           !write(99,*) -real(n-nshift)/real(nspacing),dosVB(n)/Ac

        !end do

        do n=nomega,nshiftVB,-1

           write(99,*) -real(n-nshiftVB)/real(nspacing),dosVB(n)*sqrt(3._dp)/2_dp*ndim/4_dp/Ac

        end do

        do n=nshiftCB,nomega

           write(99,*) real(n-nshiftCB)/real(nspacing),dosCB(n)*sqrt(3._dp)/2_dp*ndim/4_dp/Ac

        end do

        close(99)


end subroutine dosFullvc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dosFull(bandsTB,dos,fmu,npointsBack,ndimEV,numk,ncount,vq1,vq12,nomega,nspacing,nshift,nbounds)

 integer(dp), intent(in)    :: ndimEV,numk,ncount
 integer(dp), intent(in)    :: nomega,nspacing,nshift,nbounds(:)
 integer(dp), intent(in)    :: npointsBack(numk,numk)
 real (dp)   , intent(in)    :: fmu,vq1(2),vq12(2),bandsTB(ndimEV,ncount)
 real(dp), intent(out) :: dos(:)


 real(dp), allocatable    :: temp(:),PivotEnergy(:,:),PivotMomentum(:,:)


 integer(dp) :: n,i,ivk,ivk1,ivk2
 real(dp) :: Ac,sigma0
 real(dp) :: bandsFull(ndimEV,numk+1,numk+1)
 real(dp) :: vEnergy(ndimEV,numk+1,2),vMomentum(2,numk+1,2)




 allocate(temp(nomega))

 allocate(pivotEnergy(ndimEV,4),pivotMomentum(2,4))


 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!! Periodic continuation
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 do ivk1=0,numk-1
 do ivk2=0,numk-1
 bandsFull(:,ivk1+1,ivk2+1)=bandsTB(:,npointsBack(ivk1+1,ivk2+1))-fmu
 enddo
 enddo

 do ivk=0,numk-1
 bandsFull(:,ivk+1,numk+1)=bandsFull(:,ivk+1,1)
 bandsFull(:,numk+1,ivk+1)=bandsFull(:,1,ivk+1)
 enddo
 bandsFull(:,numk+1,numk+1)=bandsFull(:,1,1)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

dos=0._dp

do ivk1=0,numk-1
 do ivk2 = 0,numk
   vEnergy(:,ivk2+1,1)=bandsFull(:,ivk1+1,ivk2+1)
   vEnergy(:,ivk2+1,2)=bandsFull(:,ivk1+2,ivk2+1)
   vMomentum(:,ivk2+1,1)=(ivk1*vq1(:)+ivk2*vq12(:))/real(numk,dp)
   vMomentum(:,ivk2+1,2)=((ivk1+1)*vq1(:)+ivk2*vq12(:))/real(numk,dp)
 end do

 temp=0._dp

!$omp parallel do &
!$omp private(ivk2,pivotEnergy,pivotMomentum,temp) &
!$omp shared(vEnergy,vMomentum,nshift,nomega,nbounds) &
!$omp reduction(+:dos)
!!$omp reduction(+:dosJoint)
do ivk2=0,numk-1

pivotEnergy=reshape([vEnergy(:,ivk2+1,1),vEnergy(:,ivk2+1,2),vEnergy(:,ivk2+2,1),vEnergy(:,ivk2+2,2)],[ndimEV,4_dp])

pivotMomentum=reshape([vMomentum(:,ivk2+1,1),vMomentum(:,ivk2+1,2),vMomentum(:,ivk2+2,1),vMomentum(:,ivk2+2,2)],[2_dp,4_dp])

call dos_S(temp,pivotMomentum,pivotEnergy,nomega,nspacing,nshift,nbounds)

dos(:)=dos(:)+temp(:)

enddo ! ivk2
!$omp end parallel do

enddo ! ivk1



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

Ac=4._dp*pi*pi ! normalization: int DOS*sqrt(3)/2*ndim/4= number of bands =iu-il+1


        do n=1,nomega

           write(99,*) real(n-nshift)/real(nspacing),dos(n)*sqrt(3._dp)/2_dp*ndim/4_dp/Ac

        end do

        close(99)


end subroutine dosFull


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%^M

subroutine dos_S(dos,vk,energia,nomega,nspacing,nshift,nbounds)


 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: vk(:,:),energia(:,:)
 real(dp)   , intent(out) :: dos(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4)
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: add,weight


   dos=0._dp

   do i_c = nbounds(1),nbounds(4)

   energy(:)=energia(i_c,:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   call sort(energy1,energy2,energy3,vk1,vk2,vk3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift


    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         dos(iE)=dos(iE)+add
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        dos(iE)=dos(iE)+add
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   call sort(energy1,energy2,energy3,vk1,vk2,vk3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         dos(iE)=dos(iE)+add
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        dos(iE)=dos(iE)+add
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    end do

    !write(*,*) 'TestDos',maxval(dos)

end subroutine dos_S



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%^M

subroutine drude_S(drude,vk,energia,nomega,nspacing,nshift,nbounds,nxy,dk)


 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4),nxy
 real(dp)   , intent(in) :: vk(:,:),energia(:,:),dk(:)
 real(dp)   , intent(out) :: drude(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4),velocity
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: add,weight


   drude=0._dp

   do i_c = nbounds(1),nbounds(4)

   energy(:)=energia(i_c,:)

   if(nxy.EQ.1)then
    velocity=abs(energy(4)-energy(1))/dk(nxy)
   else
    velocity=abs(energy(3)-energy(2))/dk(nxy)
   endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   call sort(energy1,energy2,energy3,vk1,vk2,vk3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))*velocity

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift


    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         drude(iE)=drude(iE)+add
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        drude(iE)=drude(iE)+add
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   call sort(energy1,energy2,energy3,vk1,vk2,vk3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))*velocity

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         drude(iE)=drude(iE)+add
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        drude(iE)=drude(iE)+add
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    end do

    !write(*,*) 'TestDos',maxval(dos)

end subroutine drude_S

!!!!!!!!!!!!!!!!!!!!!!!!

subroutine drude_V(drudeX,drudeY,vk,energia,Vx,Vy,nomega,nspacing,nshift,nbounds)

 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: vk(:,:),energia(:,:),Vx(:,:),Vy(:,:)
 real(dp)   , intent(out) :: drudeX(:),drudeY(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4)
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: f1(2),f2(2),f3(2),add,weight


   drudeX=0._dp
   drudeY=0._dp

   do i_c = nbounds(1),nbounds(4)

   energy(:)=energia(i_c,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=Vx(i_c,1)
   f2(1)=Vx(i_c,2)
   f3(1)=Vx(i_c,3)

   f1(2)=Vy(i_c,1)
   f2(2)=Vy(i_c,2)
   f3(2)=Vy(i_c,3)




   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift


    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         drudeX(iE)=drudeX(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         drudeY(iE)=drudeY(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        drudeX(iE)=drudeX(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        drudeY(iE)=drudeY(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=Vx(i_c,4)
   f2(1)=Vx(i_c,2)
   f3(1)=Vx(i_c,3)

   f1(2)=Vy(i_c,4)
   f2(2)=Vy(i_c,2)
   f3(2)=Vy(i_c,3)

   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         drudeX(iE)=drudeX(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         drudeY(iE)=drudeY(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        drudeX(iE)=drudeX(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        drudeY(iE)=drudeY(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    end do

    !write(*,*) 'TestDos',maxval(dos)

end subroutine drude_V

!!!!!!!!!!!!!!!!!!!!!!!!

subroutine drude_F(drude,vk,energia,Vx,Vy,nomega,nspacing,nshift,fmu,nbounds)

 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: vk(:,:),energia(:,:),Vx(:,:),Vy(:,:),fmu
 real(dp)   , intent(out) :: drude(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4)
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: f1(1),f2(1),f3(1),add,weight


   drude=0._dp

   do i_c = nbounds(1),nbounds(4)

   energy(:)=energia(i_c,:)-fmu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=(Vx(i_c,1)*Vx(i_c,1)+Vy(i_c,1)*Vy(i_c,1))/2._dp
   f2(1)=(Vx(i_c,2)*Vx(i_c,2)+Vy(i_c,2)*Vy(i_c,2))/2._dp
   f3(1)=(Vx(i_c,2)*Vx(i_c,2)+Vy(i_c,2)*Vy(i_c,2))/2._dp





   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift


    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         drude(iE)=drude(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        drude(iE)=drude(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=(Vx(i_c,4)*Vx(i_c,4)+Vy(i_c,4)*Vy(i_c,4))/2._dp
   f2(1)=(Vx(i_c,2)*Vx(i_c,2)+Vy(i_c,2)*Vy(i_c,2))/2._dp
   f3(1)=(Vx(i_c,3)*Vx(i_c,3)+Vy(i_c,3)*Vy(i_c,3))/2._dp


   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         drude(iE)=drude(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        drude(iE)=drude(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    end do

    !write(*,*) 'TestDos',maxval(dos)

end subroutine drude_F

!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!

subroutine drudeCB(Dtot,Dmag,Dxy,vk,energia,vx1,vy1,vx2,vy2,vxI,vyI,nomega,nspacing,nshift,efermi,nbounds)

 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: efermi,vk(:,:),energia(:,:),vx1(:,:),vy1(:,:),vx2(:,:),vy2(:,:),vxI(:,:),vyI(:,:)
 real(dp)   , intent(out) :: Dtot(:),Dmag(:),Dxy(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4)
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: f1(3),f2(3),f3(3),fvx1(4),fvx2(4),fvxI(4),fvy1(4),fvy2(4),fvyI(4),add,weight


   Dtot=0._dp
   Dmag=0._dp
   Dxy=0._dp

   do i_c = nbounds(3),nbounds(4)

   energy(:)=energia(i_c,:)-efermi

   fvx1(:)=vx1(i_c,:)

   fvx2(:)=vx2(i_c,:)

   fvxI(:)=vxI(i_c,:)

   fvy1(:)=vy1(i_c,:)

   fvy2(:)=vy2(i_c,:)

   fvyI(:)=vyI(i_c,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=((fvx1(1)+fvx2(1)+fvxI(1))*(fvx1(1)+fvx2(1)+fvxI(1))+(fvy1(1)+fvy2(1)+fvyI(1))*(fvy1(1)+fvy2(1))+fvyI(1))/2._dp
   f2(1)=((fvx1(2)+fvx2(2)+fvxI(2))*(fvx1(2)+fvx2(2)+fvxI(2))+(fvy1(2)+fvy2(2)+fvyI(2))*(fvy1(2)+fvy2(2))+fvyI(2))/2._dp
   f3(1)=((fvx1(3)+fvx2(3)+fvxI(3))*(fvx1(3)+fvx2(3)+fvxI(3))+(fvy1(3)+fvy2(3)+fvyI(3))*(fvy1(3)+fvy2(3))+fvyI(3))/2._dp

   f1(2)=((fvx1(1)-fvx2(1))*(fvx1(1)-fvx2(1))+(fvy1(1)-fvy2(1))*(fvy1(1)-fvy2(1)))/2._dp
   f2(2)=((fvx1(2)-fvx2(2))*(fvx1(2)-fvx2(2))+(fvy1(2)-fvy2(2))*(fvy1(2)-fvy2(2)))/2._dp
   f3(2)=((fvx1(3)-fvx2(3))*(fvx1(3)-fvx2(3))+(fvy1(3)-fvy2(3))*(fvy1(3)-fvy2(3)))/2._dp

   f1(3)=(fvx1(1)*(fvy2(1)+fvyI(1))-(fvx2(1)+fvxI(1))*fvy1(1))/2._dp
   f2(3)=(fvx1(2)*(fvy2(2)+fvyI(2))-(fvx2(2)+fvxI(2))*fvy1(2))/2._dp
   f3(3)=(fvx1(3)*(fvy2(3)+fvyI(3))-(fvx2(3)+fvxI(3))*fvy1(3))/2._dp




   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift


    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         Dtot(iE)=Dtot(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         Dmag(iE)=Dmag(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
         Dxy(iE)=Dxy(iE)+add*(f3(3)*(n2-iE)/real(n2-n3)+f2(3)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        Dtot(iE)=Dtot(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        Dmag(iE)=Dmag(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
        Dxy(iE)=Dxy(iE)+add*(f2(3)*(n1-iE)/real(n1-n2)+f1(3)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=((fvx1(4)+fvx2(4)+fvxI(4))*(fvx1(4)+fvx2(4)+fvxI(4))+(fvy1(4)+fvy2(4)+fvyI(4))*(fvy1(4)+fvy2(4))+fvyI(4))/2._dp
   f2(1)=((fvx1(2)+fvx2(2)+fvxI(2))*(fvx1(2)+fvx2(2)+fvxI(2))+(fvy1(2)+fvy2(2)+fvyI(2))*(fvy1(2)+fvy2(2))+fvyI(2))/2._dp
   f3(1)=((fvx1(3)+fvx2(3)+fvxI(3))*(fvx1(3)+fvx2(3)+fvxI(3))+(fvy1(3)+fvy2(3)+fvyI(3))*(fvy1(3)+fvy2(3))+fvyI(3))/2._dp

   f1(2)=((fvx1(4)-fvx2(4))*(fvx1(4)-fvx2(4))+(fvy1(4)-fvy2(4))*(fvy1(4)-fvy2(4)))/2._dp
   f2(2)=((fvx1(2)-fvx2(2))*(fvx1(2)-fvx2(2))+(fvy1(2)-fvy2(2))*(fvy1(2)-fvy2(2)))/2._dp
   f3(2)=((fvx1(3)-fvx2(3))*(fvx1(3)-fvx2(3))+(fvy1(3)-fvy2(3))*(fvy1(3)-fvy2(3)))/2._dp

   f1(3)=(fvx1(4)*(fvy2(4)+fvyI(4))-(fvx2(4)+fvxI(4))*fvy1(4))/2._dp
   f2(3)=(fvx1(2)*(fvy2(2)+fvyI(2))-(fvx2(2)+fvxI(2))*fvy1(2))/2._dp
   f3(3)=(fvx1(3)*(fvy2(3)+fvyI(3))-(fvx2(3)+fvxI(3))*fvy1(3))/2._dp

   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         Dtot(iE)=Dtot(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         Dmag(iE)=Dmag(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
         Dxy(iE)=Dxy(iE)+add*(f3(3)*(n2-iE)/real(n2-n3)+f2(3)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        Dtot(iE)=Dtot(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        Dmag(iE)=Dmag(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
        Dxy(iE)=Dxy(iE)+add*(f2(3)*(n1-iE)/real(n1-n2)+f1(3)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    end do

    !write(*,*) 'TestDos',maxval(dos)

end subroutine drudeCB

!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!

subroutine drudeVB(Dtot,Dmag,Dxy,vk,energia,vx1,vy1,vx2,vy2,vxI,vyI,nomega,nspacing,nshift,efermi,nbounds)

 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: efermi,vk(:,:),energia(:,:),vx1(:,:),vy1(:,:),vx2(:,:),vy2(:,:),vxI(:,:),vyI(:,:)
 real(dp)   , intent(out) :: Dtot(:),Dmag(:),Dxy(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4)
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: f1(3),f2(3),f3(3),fvx1(4),fvx2(4),fvxI(4),fvy1(4),fvy2(4),fvyI(4),add,weight


   Dtot=0._dp
   Dmag=0._dp
   Dxy=0._dp

   do i_v = nbounds(1),nbounds(2)

   energy(:)=efermi-energia(i_v,:)

   fvx1(:)=vx1(i_v,:)

   fvx2(:)=vx2(i_v,:)

   fvxI(:)=vxI(i_v,:)

   fvy1(:)=vy1(i_v,:)

   fvy2(:)=vy2(i_v,:)

   fvyI(:)=vyI(i_v,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=((fvx1(1)+fvx2(1)+fvxI(1))*(fvx1(1)+fvx2(1)+fvxI(1))+(fvy1(1)+fvy2(1)+fvyI(1))*(fvy1(1)+fvy2(1))+fvyI(1))/2._dp
   f2(1)=((fvx1(2)+fvx2(2)+fvxI(2))*(fvx1(2)+fvx2(2)+fvxI(2))+(fvy1(2)+fvy2(2)+fvyI(2))*(fvy1(2)+fvy2(2))+fvyI(2))/2._dp
   f3(1)=((fvx1(3)+fvx2(3)+fvxI(3))*(fvx1(3)+fvx2(3)+fvxI(3))+(fvy1(3)+fvy2(3)+fvyI(3))*(fvy1(3)+fvy2(3))+fvyI(3))/2._dp

   f1(2)=((fvx1(1)-fvx2(1))*(fvx1(1)-fvx2(1))+(fvy1(1)-fvy2(1))*(fvy1(1)-fvy2(1)))/2._dp
   f2(2)=((fvx1(2)-fvx2(2))*(fvx1(2)-fvx2(2))+(fvy1(2)-fvy2(2))*(fvy1(2)-fvy2(2)))/2._dp
   f3(2)=((fvx1(3)-fvx2(3))*(fvx1(3)-fvx2(3))+(fvy1(3)-fvy2(3))*(fvy1(3)-fvy2(3)))/2._dp

   f1(3)=(fvx1(1)*(fvy2(1)+fvyI(1))-(fvx2(1)+fvxI(1))*fvy1(1))/2._dp
   f2(3)=(fvx1(2)*(fvy2(2)+fvyI(2))-(fvx2(2)+fvxI(2))*fvy1(2))/2._dp
   f3(3)=(fvx1(3)*(fvy2(3)+fvyI(3))-(fvx2(3)+fvxI(3))*fvy1(3))/2._dp




   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift


    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         Dtot(iE)=Dtot(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         Dmag(iE)=Dmag(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
         Dxy(iE)=Dxy(iE)+add*(f3(3)*(n2-iE)/real(n2-n3)+f2(3)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        Dtot(iE)=Dtot(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        Dmag(iE)=Dmag(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
        Dxy(iE)=Dxy(iE)+add*(f2(3)*(n1-iE)/real(n1-n2)+f1(3)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=((fvx1(4)+fvx2(4)+fvxI(4))*(fvx1(4)+fvx2(4)+fvxI(4))+(fvy1(4)+fvy2(4)+fvyI(4))*(fvy1(4)+fvy2(4))+fvyI(4))/2._dp
   f2(1)=((fvx1(2)+fvx2(2)+fvxI(2))*(fvx1(2)+fvx2(2)+fvxI(2))+(fvy1(2)+fvy2(2)+fvyI(2))*(fvy1(2)+fvy2(2))+fvyI(2))/2._dp
   f3(1)=((fvx1(3)+fvx2(3)+fvxI(3))*(fvx1(3)+fvx2(3)+fvxI(3))+(fvy1(3)+fvy2(3)+fvyI(3))*(fvy1(3)+fvy2(3))+fvyI(3))/2._dp

   f1(2)=((fvx1(4)-fvx2(4))*(fvx1(4)-fvx2(4))+(fvy1(4)-fvy2(4))*(fvy1(4)-fvy2(4)))/2._dp
   f2(2)=((fvx1(2)-fvx2(2))*(fvx1(2)-fvx2(2))+(fvy1(2)-fvy2(2))*(fvy1(2)-fvy2(2)))/2._dp
   f3(2)=((fvx1(3)-fvx2(3))*(fvx1(3)-fvx2(3))+(fvy1(3)-fvy2(3))*(fvy1(3)-fvy2(3)))/2._dp

   f1(3)=(fvx1(4)*(fvy2(4)+fvyI(4))-(fvx2(4)+fvxI(4))*fvy1(4))/2._dp
   f2(3)=(fvx1(2)*(fvy2(2)+fvyI(2))-(fvx2(2)+fvxI(2))*fvy1(2))/2._dp
   f3(3)=(fvx1(3)*(fvy2(3)+fvyI(3))-(fvx2(3)+fvxI(3))*fvy1(3))/2._dp

   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         Dtot(iE)=Dtot(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         Dmag(iE)=Dmag(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
         Dxy(iE)=Dxy(iE)+add*(f3(3)*(n2-iE)/real(n2-n3)+f2(3)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        Dtot(iE)=Dtot(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        Dmag(iE)=Dmag(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
        Dxy(iE)=Dxy(iE)+add*(f2(3)*(n1-iE)/real(n1-n2)+f1(3)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    end do

    !write(*,*) 'TestDos',maxval(dos)

end subroutine drudeVB

!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!

subroutine drudeFull(Dtot,Dmag,Dxy,vk,energia,vx1,vy1,vx2,vy2,vxI,vyI,nomega,nspacing,nshift,efermi,nbounds)

 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: efermi,vk(:,:),energia(:,:),vx1(:,:),vy1(:,:),vx2(:,:),vy2(:,:),vxI(:,:),vyI(:,:)
 real(dp)   , intent(out) :: Dtot(:),Dmag(:),Dxy(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4)
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: f1(3),f2(3),f3(3),fvx1(4),fvx2(4),fvxI(4),fvy1(4),fvy2(4),fvyI(4),add,weight


   Dtot=0._dp
   Dmag=0._dp
   Dxy=0._dp

   do i_c = nbounds(1),nbounds(4)

   energy(:)=energia(i_c,:)-efermi

   fvx1(:)=vx1(i_c,:)

   fvx2(:)=vx2(i_c,:)

   fvxI(:)=vxI(i_c,:)

   fvy1(:)=vy1(i_c,:)

   fvy2(:)=vy2(i_c,:)

   fvyI(:)=vyI(i_c,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=((fvx1(1)+fvx2(1)+fvxI(1))*(fvx1(1)+fvx2(1)+fvxI(1))+(fvy1(1)+fvy2(1)+fvyI(1))*(fvy1(1)+fvy2(1))+fvyI(1))/2._dp
   f2(1)=((fvx1(2)+fvx2(2)+fvxI(2))*(fvx1(2)+fvx2(2)+fvxI(2))+(fvy1(2)+fvy2(2)+fvyI(2))*(fvy1(2)+fvy2(2))+fvyI(2))/2._dp
   f3(1)=((fvx1(3)+fvx2(3)+fvxI(3))*(fvx1(3)+fvx2(3)+fvxI(3))+(fvy1(3)+fvy2(3)+fvyI(3))*(fvy1(3)+fvy2(3))+fvyI(3))/2._dp

   f1(2)=((fvx1(1)-fvx2(1))*(fvx1(1)-fvx2(1))+(fvy1(1)-fvy2(1))*(fvy1(1)-fvy2(1)))/2._dp
   f2(2)=((fvx1(2)-fvx2(2))*(fvx1(2)-fvx2(2))+(fvy1(2)-fvy2(2))*(fvy1(2)-fvy2(2)))/2._dp
   f3(2)=((fvx1(3)-fvx2(3))*(fvx1(3)-fvx2(3))+(fvy1(3)-fvy2(3))*(fvy1(3)-fvy2(3)))/2._dp

   f1(3)=(fvx1(1)*(fvy2(1)+fvyI(1))-(fvx2(1)+fvxI(1))*fvy1(1))/2._dp
   f2(3)=(fvx1(2)*(fvy2(2)+fvyI(2))-(fvx2(2)+fvxI(2))*fvy1(2))/2._dp
   f3(3)=(fvx1(3)*(fvy2(3)+fvyI(3))-(fvx2(3)+fvxI(3))*fvy1(3))/2._dp




   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift


    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         Dtot(iE)=Dtot(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         Dmag(iE)=Dmag(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
         Dxy(iE)=Dxy(iE)+add*(f3(3)*(n2-iE)/real(n2-n3)+f2(3)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        Dtot(iE)=Dtot(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        Dmag(iE)=Dmag(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
        Dxy(iE)=Dxy(iE)+add*(f2(3)*(n1-iE)/real(n1-n2)+f1(3)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=((fvx1(4)+fvx2(4)+fvxI(4))*(fvx1(4)+fvx2(4)+fvxI(4))+(fvy1(4)+fvy2(4)+fvyI(4))*(fvy1(4)+fvy2(4))+fvyI(4))/2._dp
   f2(1)=((fvx1(2)+fvx2(2)+fvxI(2))*(fvx1(2)+fvx2(2)+fvxI(2))+(fvy1(2)+fvy2(2)+fvyI(2))*(fvy1(2)+fvy2(2))+fvyI(2))/2._dp
   f3(1)=((fvx1(3)+fvx2(3)+fvxI(3))*(fvx1(3)+fvx2(3)+fvxI(3))+(fvy1(3)+fvy2(3)+fvyI(3))*(fvy1(3)+fvy2(3))+fvyI(3))/2._dp

   f1(2)=((fvx1(4)-fvx2(4))*(fvx1(4)-fvx2(4))+(fvy1(4)-fvy2(4))*(fvy1(4)-fvy2(4)))/2._dp
   f2(2)=((fvx1(2)-fvx2(2))*(fvx1(2)-fvx2(2))+(fvy1(2)-fvy2(2))*(fvy1(2)-fvy2(2)))/2._dp
   f3(2)=((fvx1(3)-fvx2(3))*(fvx1(3)-fvx2(3))+(fvy1(3)-fvy2(3))*(fvy1(3)-fvy2(3)))/2._dp

   f1(3)=(fvx1(4)*(fvy2(4)+fvyI(4))-(fvx2(4)+fvxI(4))*fvy1(4))/2._dp
   f2(3)=(fvx1(2)*(fvy2(2)+fvyI(2))-(fvx2(2)+fvxI(2))*fvy1(2))/2._dp
   f3(3)=(fvx1(3)*(fvy2(3)+fvyI(3))-(fvx2(3)+fvxI(3))*fvy1(3))/2._dp

   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         Dtot(iE)=Dtot(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         Dmag(iE)=Dmag(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
         Dxy(iE)=Dxy(iE)+add*(f3(3)*(n2-iE)/real(n2-n3)+f2(3)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        Dtot(iE)=Dtot(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        Dmag(iE)=Dmag(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
        Dxy(iE)=Dxy(iE)+add*(f2(3)*(n1-iE)/real(n1-n2)+f1(3)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    end do

    !write(*,*) 'TestDos',maxval(dos)

end subroutine drudeFull



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dosVB_S(dos,vk,energia,nomega,nspacing,nshift,efermi,nbounds)


 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: efermi
 real(dp)   , intent(in) :: vk(:,:),energia(:,:)
 real(dp)   , intent(out) :: dos(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE
 real(dp) :: energy(4)
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: add,weight

 dos=0._dp

 do i_v = nbounds(1),nbounds(2)

   energy(:)=efermi-energia(i_v,:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   call sort(energy1,energy2,energy3,vk1,vk2,vk3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         dos(iE)=dos(iE)+add
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        dos(iE)=dos(iE)+add
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   call sort(energy1,energy2,energy3,vk1,vk2,vk3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         dos(iE)=dos(iE)+add
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        dos(iE)=dos(iE)+add
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    end do


end subroutine dosVB_S

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dosCB_S(dos,vk,energia,nomega,nspacing,nshift,efermi,nbounds)


 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: efermi
 real(dp)   , intent(in) :: vk(:,:),energia(:,:)
 real(dp)   , intent(out) :: dos(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE
 real(dp) :: energy(4)
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: add,weight

 dos=0._dp

 do i_v = nbounds(3),nbounds(4)

   energy(:)=energia(i_v,:)-efermi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   call sort(energy1,energy2,energy3,vk1,vk2,vk3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         dos(iE)=dos(iE)+add
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        dos(iE)=dos(iE)+add
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   call sort(energy1,energy2,energy3,vk1,vk2,vk3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         dos(iE)=dos(iE)+add
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        dos(iE)=dos(iE)+add
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    end do


end subroutine dosCB_S

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dosJoint_S(dos,vk,energia,nomega,nspacing,nshift,nbounds)


 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: vk(:,:),energia(:,:)
 real(dp)   , intent(out) :: dos(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE
 real(dp) :: energy(4)
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: add,weight

 dos=0._dp

 do i_v = nbounds(1),nbounds(2)
   do i_c = nbounds(3),nbounds(4)


   energy(:)=energia(i_c,:)-energia(i_v,:)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   call sort(energy1,energy2,energy3,vk1,vk2,vk3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         dos(iE)=dos(iE)+add
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        dos(iE)=dos(iE)+add
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   call sort(energy1,energy2,energy3,vk1,vk2,vk3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         dos(iE)=dos(iE)+add
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        dos(iE)=dos(iE)+add
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    end do
    end do

end subroutine dosJoint_S

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



subroutine ComplexH0_delta_magneticfield(zh,coord,ndimq,vk,Delta,tn,vq1,vq2,p,q)

  integer(dp), intent(in)    :: ndimq, p, q
  real(dp)   , intent(in)    :: tn(2,3), vk(2), vq1(2), vq2(2), coord(ndimq,3),Delta
  complex(dp), intent(inout) :: zh(ndimq,ndimq)

  integer(dp) :: i, j, nt
  real(dp)    :: rij, coord_diff
  complex(dp) :: zi_vk_tn, ph

  zh(:,:) = cmplx(0.0_dp,0.0_dp,dp)

  ! We are always doing Delta=0, so this doesn't make a difference
  do i = 1 , ndimq/q/4
    do j = 0,q-1
      zh(i+j*ndimq/q,i+j*ndimq/q) = zh(i+j*ndimq/q,i+j*ndimq/q) + Delta !A1
      !zh(i+ndimq/q/2+j*ndimq/q,i+ndimq/q/2+j*ndimq/q) = zh(i+ndimq/q/2+j*ndimq/q,i+ndimq/q/2+j*ndimq/q) + Delta  !A2
      zh(i+ndimq/q/4+j*ndimq/q,i+ndimq/q/4+j*ndimq/q) = zh(i+ndimq/q/4+j*ndimq/q,i+ndimq/q/4+j*ndimq/q) - Delta !B1
      !zh(i+3*ndimq/q/4+j*ndimq/q,i+3*ndimq/q/4+j*ndimq/q) = zh(i+3*ndimq/q/4+j*ndimq/q,i+3*ndimq/q/4+j*ndimq/q) - Delta  !B2
    enddo
  enddo

  do i = 1 , ndimq
    do j = 1 , ndimq

      ! Compute the phase factor(we have passed the magnetic lattice vectors(tn) and the Moire reciprocal vectors(vq).
      ! Also, notice that in the last argument ntheta=0.0_dp.)
      rij = norm2(coord(i,:)-coord(j,:))
      coord_diff = (coord(i,3)-coord(j,3))**2
      ph = phase([coord(i,1:2),coord(j,1:2)],[vq1, vq2],[tn(1:2,1), tn(1:2,2)],p,q)
      ! ph = phase2([coord(i,1:2),coord(j,1:2)],[vq1, vq2],[tn(1:2,1), tn(1:2,2)],p,q)


      if(coord_diff < deltaR) then
        zh(i,j) = zh(i,j) + 0.5_dp*ft(rij)*ph
        zh(j,i) = zh(j,i) + 0.5_dp*ft(rij)*conjg(ph)
      else
        zh(i,j) = zh(i,j) + 0.5_dp*ftperp(rij)*ph
        zh(j,i) = zh(j,i) + 0.5_dp*ftperp(rij)*conjg(ph)
      endif

      !!!!!! Coupling to adjacent Wigner-Seitz cell with tn
      do nt = 1 , 3
        rij = norm2(coord(i,:)-[coord(j,1:2)-tn(:,nt),coord(j,3)])
        zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
        ph = phase([coord(i,1:2),coord(j,1:2)-tn(:,nt)],[vq1, vq2],[tn(:,1), tn(:,2)],p,q)
        ! ph = phase2([coord(i,1:2),coord(j,1:2)-tn(:,nt)],[vq1, vq2],[tn(:,1), tn(:,2)],p,q)

        if(coord_diff < deltaR)then
          zh(i,j)=zh(i,j)+ft(rij)*exp( zi_vk_tn)*ph
          ! zh(j,i)=zh(j,i)+ft(rij)*exp(-zi_vk_tn)*phase([coord(j,1:2)-tn(:,nt),coord(i,1:2)],[vq1, vq2],[tn(:,1), tn(:,2)],p,q)

          zh(j,i)=zh(j,i) + ft(rij)*exp(-zi_vk_tn)*conjg(ph)

        else
          zh(i,j)=zh(i,j)+ftperp(rij)*exp( zi_vk_tn)*ph
          ! zh(j,i)=zh(j,i)+ftperp(rij)*exp(-zi_vk_tn)* &
            ! phase([coord(j,1:2)-tn(:,nt),coord(i,1:2)],[vq1, vq2],[tn(:,1), tn(:,2)],p,q)

          zh(j,i)=zh(j,i)+ ftperp(rij)*exp(-zi_vk_tn)*conjg(ph)
        endif

      end do

    end do
  end do

 end subroutine ComplexH0_delta_magneticfield

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function phase(coord,vq,ts,p,q)
  integer(dp), intent(in)    :: p, q
  real(dp)   , intent(in)    :: ts(4), coord(4), vq(4)
  ! complex(dp), intent(inout) :: zh

  complex(dp) :: phase
  integer(dp) :: i, j
  real(dp)    :: ei(2), ej(2), ci(2), cj(2), vq1(2), vq2(2)
  complex(dp) :: zi = (0.0_dp,1.0_dp), ph01, ph02, ph0

  vq1 = vq(1:2)
  vq2 = vq(3:4)
  ci = coord(1:2)
  cj = coord(3:4)
  ! write(*,*) ci
  ! write(*,*) cj

  ! Coordinates in terms of t1 and t2(the Moire lattice vectors. Notice that we passed tq, so we have t2=t2q/q)
  call ecoords(ei,ci,[ts(1:2), ts(3:4)/real(q)])
  call ecoords(ej,cj,[ts(1:2), ts(3:4)/real(q)])

  ! Compute phases. See the definitions of rho and sumek below
  if (abs(ej(1)-ei(1)) .gt. 1e-9_dp) then
    ph01 = dot_product(vq2,cj-ci)/(ej(1)-ei(1))* &
    ((ej(1)**2-ei(1)**2)/2 - (rho(ej(1)) - rho(ei(1))))
    ph02 = dot_product(vq1,cj-ci)/abs(ej(1)-ei(1))*sumek(ei,ej)
    ph0 = ph01 - ph02

    phase = exp(zi*ph0*real(p,dp)/real(q,dp))

  else
    ph0 = (ei(1) - real(floor(ei(1)+1e-9_dp),dp))*dot_product(vq2,cj-ci)
    phase = exp(zi*ph0*real(p,dp)/real(q,dp))
  endif

  end function phase

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function phase2(coord,vq,ts,p,q)

    integer(dp), intent(in)    :: p, q
    real(dp)   , intent(in)    :: ts(4), coord(4), vq(4)
    ! complex(dp), intent(inout) :: zh

    complex(dp) :: phase2
    integer(dp) :: i, j
    real(dp)    :: ei(2), ej(2), ci(2), cj(2), vq1(2), vq2(2)
    complex(dp) :: zi = (0.0_dp,1.0_dp), ph0

    vq1 = vq(1:2)
    vq2 = vq(3:4)
    ci = coord(1:2)
    cj = coord(3:4)
    ! write(*,*) ci
    ! write(*,*) cj

    ! Coordinates in terms of t1 and t2(the Moire lattice vectors. Notice that we passed tq, so we have t2=t2q/q)
    call ecoords(ei,ci,[ts(1:2), ts(3:4)/real(q)])
    call ecoords(ej,cj,[ts(1:2), ts(3:4)/real(q)])

    ! Compute phases. See the definitions of rho and sumek below
    ph0 = dot_product(vq2,cj-ci)*(ej(1) + ei(1))/2 - 2*pi*ej(2)*floor(ej(1)+1e-9) +  2*pi*ei(2)*floor(ei(1)+1e-9)
    phase2 = exp(zi*ph0*real(p,dp)/real(q,dp))


  end function phase2

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine ecoords(e,c,ts)
    real(dp), intent(in) :: c(2), ts(4)
    real(dp), intent(inout) :: e(2)
    ! Calculate e such that e1*t1 + e2*t2 = c
    ! Solve linear equation (ts(1)=t1(1), ts(2)=t1(2), ts(3)=t2(1), ts(4)=t2(2))
    e = 1/(ts(1)*ts(4)-ts(2)*ts(3))*[ts(4)*c(1)-ts(3)*c(2), -ts(2)*c(1)+ts(1)*c(2)]
  end subroutine ecoords

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function sumek(ei,ej)
    real(dp), intent(in) :: ei(2), ej(2)

    real(dp) :: sumek
    integer(4) :: i, d
    ! Sum of the different values of e2 when e1=integer along the path
    sumek = 0.0_dp
    ! d tells us the number of times e1=integer along the path (notice +1e-9 to avoid rounding errors)
    d = floor(ej(1)+1e-9_dp) - floor(ei(1)+1e-9_dp)
    if (d .ne. 0) then
      ! if d>0, we loop from floor(ei(1))+1 to floor(ej(1)). Else we loop from floor(ei(1)) to floor(ej(1))+1
      do i=floor(ei(1)+1e-9_dp)+(d/iabs(d)+1)/2, floor(ei(1)+1e-9_dp)+d-(d/iabs(d)-1)/2, d/iabs(d)
        !we add the value of e2 when e1=i
        sumek = sumek + ei(2) + (ei(1)-real(i,dp))*(ej(2)-ei(2))/(ei(1)-ej(1))
      enddo
    endif
  end function sumek

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function rho(x)
    real(dp), intent(in) :: x
    real(dp) :: rho
    ! Notice the -1e-9 to avoid rounding errors when x is exactly an integer
    rho = real(floor(x+1e-9_dp),dp)*real(floor(x+1e-9_dp)-1,dp)/2 + real(floor(x+1e-9_dp),dp)*(x - real(floor(x+1e-9_dp),dp))

  end function rho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function gcd(a,b)
    ! We compute the greatest common divisor by Euclid's method
    integer(dp), intent(in) :: a,b

    integer(dp) :: a1,b1,c,gcd
    a1 = a
    b1 = b
    if ((a1.EQ.0).OR.(b1.EQ.0)) then
      if ((a1.EQ.1).OR.(b1.EQ.1)) then
        b1 = 1
      else
        b1 = 0
      endif
    else
      if (b1 .gt. a1) then
        c=a1
        a1=b1
        b1=c
      endif
      do
        c = mod(a1,b1)
        if (c .eq. 0) exit
        a1 = b1
        b1 = c
      enddo
    endif
    gcd = b1

  end function gcd

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine vx_magneticfield(zvx,ndim,coord,vk,tn,xy,vq1,vq2,p,q)

    integer(dp), intent(in)    :: ndim, xy, p, q
    real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3), vq1(2), vq2(2)
    complex(dp), intent(inout) :: zvx(:,:)

    integer(dp) :: i, j, t
    real(dp) :: rij, rxy, coord_diff
    complex(dp) :: zi_vk_tn, ph

    zvx = cmplx(0.0_dp,0.0_dp,dp)

    do i = 1 , ndim
      do j = 1 , ndim

        rij = norm2(coord(i,:)-coord(j,:))
        rxy  = coord(i,xy)-coord(j,xy)
        coord_diff = (coord(i,3)-coord(j,3))**2
        ph = phase2([coord(i,1:2),coord(j,1:2)],[vq1, vq2],[tn(1:2,1), tn(1:2,2)],p,q)

        if(coord_diff < deltaR .and. rij > deltaR) then
          zvx(i,j) = cmplx(0.0_dp,ft(rij)*rxy,dp)*ph
        end if

        if(coord_diff > deltaR)then
          zvx(i,j) = cmplx(0.0_dp,ftperp(rij)*rxy,dp)*ph
        endif

   !!!!!!! Coupling to adjacent Wigner-Seitz cell with tn
        do t = 1 , 3
          rij = norm2(coord(i,:)-[coord(j,1:2)-tn(:,t),coord(j,3)])
          rxy=(coord(i,xy)-coord(j,xy)+tn(xy,t))
          zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,t)),dp)
          ph = phase2([coord(i,1:2),coord(j,1:2)-tn(:,t),coord(j,3)],[vq1, vq2],[tn(1:2,1), tn(1:2,2)],p,q)

          if(coord_diff < deltaR) then
            zvx(i,j)=zvx(i,j)+ft(rij)*exp( zi_vk_tn)*cmplx(0.0_dp,rxy,dp)*ph
            zvx(j,i)=zvx(j,i)-ft(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)*conjg(ph)
          else
            zvx(i,j)=zvx(i,j)+ftperp(rij)*exp( zi_vk_tn)*cmplx(0.0_dp,rxy,dp)*ph
            zvx(j,i)=zvx(j,i)-ftperp(rij)*exp(-zi_vk_tn)*cmplx(0.0_dp,rxy,dp)*conjg(ph)
          endif
        end do

      end do
    end do

   end subroutine vx_magneticfield


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine sort(energyA,energyB,energyC,vkA,vkB,vkC)

 real(dp)   , intent(inout) :: energyA,energyB,energyC
 real(dp)   , intent(inout) :: vkA(:),vkB(:),vkC(:)
 real(dp)   :: energy1,energy2,energy3
 real(dp)   :: vk1(2),vk2(2),vk3(2)


       if(energyA.GT.energyB)then
    if(energyA.GT.energyC)then
        if(energyB.GT.energyC)then
            energy1=energyA
                          energy2=energyB
                          energy3=energyC

            vk1=vkA
            vk2=vkB
                          vk3=vkC
                 else
            energy1=energyA
                          energy2=energyC
                          energy3=energyB

                          vk1=vkA
             vk2=vkC
                          vk3=vkB
        endif
    else
        energy1=energyC
                 energy2=energyA
                 energy3=energyB

                 vk1=vkC
             vk2=vkA
                 vk3=vkB
    endif
    else
    if(energyB.GT.energyC)then
        if(energyC.GT.energyA)then
            energy1=energyB
                          energy2=energyC
                          energy3=energyA

                          vk1=vkB
            vk2=vkC
                          vk3=vkA
                 else
            energy1=energyB
                          energy2=energyA
                          energy3=energyC

                          vk1=vkB
            vk2=vkA
                          vk3=vkC
        endif
    else
        energy1=energyC
                energy2=energyB
                energy3=energyA

                vk1=vkC
        vk2=vkB
                vk3=vkA
    endif
    endif

    energyA=energy1
    energyB=energy2
    energyC=energy3

    vkA=vk1
    vkB=vk2
    vkC=vk3

end subroutine sort

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine sortF(energyA,energyB,energyC,vkA,vkB,vkC,fA,fB,fC)

 real(dp)   , intent(inout) :: energyA,energyB,energyC
 real(dp)   , intent(inout) :: vkA(:),vkB(:),vkC(:)
 real(dp)   , intent(inout) :: fA(:),fB(:),fC(:)
 real(dp)   :: energy1,energy2,energy3
 real(dp)   :: vk1(2),vk2(2),vk3(2)
 real(dp) , allocatable  :: f1(:),f2(:),f3(:)

 allocate(f1(size(fA)),f2(size(fA)),f3(size(fA)))


       if(energyA.GT.energyB)then
    if(energyA.GT.energyC)then
        if(energyB.GT.energyC)then
            energy1=energyA
                        energy2=energyB
                        energy3=energyC

            vk1=vkA
            vk2=vkB
                        vk3=vkC

            f1=fA
                        f2=fB
                        f3=fC
                 else
            energy1=energyA
                          energy2=energyC
                          energy3=energyB

                          vk1=vkA
             vk2=vkC
                          vk3=vkB

                        f1=fA
                        f2=fC
                        f3=fB
        endif
    else
        energy1=energyC
                 energy2=energyA
                 energy3=energyB

                 vk1=vkC
             vk2=vkA
                 vk3=vkB

                 f1=fC
                 f2=fA
                 f3=fB
    endif
    else
    if(energyB.GT.energyC)then
        if(energyC.GT.energyA)then
            energy1=energyB
                          energy2=energyC
                          energy3=energyA

                          vk1=vkB
            vk2=vkC
                          vk3=vkA

                         f1=fB
                        f2=fC
                        f3=fA
                 else
            energy1=energyB
                          energy2=energyA
                          energy3=energyC

                          vk1=vkB
            vk2=vkA
                          vk3=vkC

                        f1=fB
                        f2=fA
                        f3=fC
        endif
    else
        energy1=energyC
                 energy2=energyB
                 energy3=energyA

                 vk1=vkC
        vk2=vkB
                 vk3=vkA

                f1=fC
                f2=fB
                f3=fA
    endif
    endif

    energyA=energy1
    energyB=energy2
    energyC=energy3

    vkA=vk1
    vkB=vk2
    vkC=vk3

    fA=f1
    fB=f2
    fC=f3

end subroutine sortF

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

subroutine overlap(rout,zev,zvxy,coord,vk,tn,ndimEV)

 integer(dp), intent(in) :: ndimEV
 real(dp)   , intent(in) :: vk(:),coord(:,:),tn(:,:)
 complex(dp), intent(in) :: zev(:,:),zvxy(:,:)
 real(dp), intent(out) :: rout(:)

 integer(dp) :: i,n






 !ztemp=0._dp

 !write(*,*) 'TestZEV',size(zev(:,1)),size(zev(1,:))

 !call vxy_MF(zvxy,ndimEV,vk,i)



 !write(*,*) 'TestTest',size(zvxy(:,1)),size(zev(:,1)),size(ztemp(1,:)),size(zout(:,1,1))

 !call matrix_overlap(zev,zvxy,ztemp)

 call z_matrix_elements_diag_real(zvxy,zev,rout)

 !call matrix_elements(zvxy,transpose(zev),rout(:,i))

 !write(*,*) 'TestZTEMP',size(ztemp(:,1)),size(ztemp(1,:))

 !call matrix_overlap(ztemp,zev,zout(:,:,i))




 !write(*,*) 'Test',rout(1,1),rout(2,1)
 !write(*,*) 'Test',rout(3,1),rout(4,1)

end subroutine overlap

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine overlapZ(zout,zev,vk,ndim)

 integer(dp), intent(in) :: ndim
 real(dp)   , intent(in) :: vk(:)
 complex(dp), intent(in) :: zev(:,:)
 complex(dp), intent(inout) :: zout(:,:,:)

 integer(dp) :: i,n
 complex(dp), allocatable :: zvxy(:,:),ztemp(:,:)

 allocate(zvxy(ndim,ndim),ztemp(ndim,ndim))


 do i=1,2

 ztemp=0._dp

 !call vxy(zvxy,ndim,vk,i)

 call matrix_overlap(zvxy,zev,ztemp)

 call matrix_overlap(zev,ztemp,zout(:,:,i))

 enddo



end subroutine overlapZ

subroutine DrudeT(dosVBin,dosCBin,drudeTotVBin,drudeTotCBin,drudeMagVBin,drudeMagCBin,drudeChiVBin,drudeChiCBin,nshift0,nomega0,nspacing,nNum)

 integer(dp), intent(in) :: nshift0,nomega0,nspacing,nNum
 real(dp)   , intent(in) :: dosVBin(:),dosCBin(:),drudeTotVBin(:),drudeTotCBin(:),drudeMagVBin(:),drudeMagCBin(:),drudeChiVBin(:),drudeChiCBin(:)

 integer(dp) :: i,n,imu,iT,iNum,nshift,nomega
 real(dp) :: A_0,A_n,Ac,vfermi,T,fmu,ome,units,THz,Lorentzian

 integer(dp) , parameter :: nstep=1

 real(dp), parameter :: deg = 2._dp
 real(dp), parameter :: degC = 4._dp
 real(dp), parameter :: t0 = 2.7_dp
 real(dp), allocatable :: occVB(:,:),occCB(:,:),occ_g(:,:)
 real(dp), allocatable :: dosVB0(:),dosCB0(:),drudeVB0(:,:),drudeCB0(:,:)
 real(dp), allocatable :: dosVB(:,:),dosCB(:,:),drudeVB(:,:,:),drudeCB(:,:,:),dos_g(:,:),drude_g(:,:)


 nshift=nshift0*12_dp
 nomega=(nomega0-nshift0)*2_dp+nshift

 allocate(occVB(nomega,nT),occCB(nomega,nT),occ_g(nomega,nT))
 allocate(dosVB0(nomega),dosCB0(nomega),drudeVB0(nomega,nNum),drudeCB0(nomega,nNum))
 allocate(dosVB(nomega,nT),dosCB(nomega,nT),drudeVB(nomega,nT,nNum),drudeCB(nomega,nT,nNum),dos_g(nomega,nT),drude_g(nomega,nT))

A_0=sqrt(3._dp)/2._dp*ndim/4._dp ! unit of density of one moire cell

A_n=2.46*2.46/10000._dp ! 10^-10/10000 = 10^-14  = 10^-12 (10^-2 m)= 10^-12 cm^-2

Ac=4._dp*pi*pi ! normalization

dosVB0=0._dp
dosCB0=0._dp

drudeVB0=0._dp
drudeCB0=0._dp


occCB=0._dp
occVB=0._dp
occ_g=0._dp


drudeVB=0._dp
drudeCB=0._dp


do n=1,nomega0-nshift0
dosVB0(n+nshift)=dosVBin(n+nshift0)
dosCB0(n+nshift)=dosCBin(n+nshift0)
drudeVB0(n+nshift,1)=drudeTotVBin(n+nshift0)
drudeCB0(n+nshift,1)=drudeTotCBin(n+nshift0)
if(nNum.GE.2)then
drudeVB0(n+nshift,2)=drudeMagVBin(n+nshift0)
drudeCB0(n+nshift,2)=drudeMagCBin(n+nshift0)
endif
if(nNum.GE.3)then
drudeVB0(n+nshift,3)=drudeChiVBin(n+nshift0)
drudeCB0(n+nshift,3)=drudeChiCBin(n+nshift0)
endif
enddo

do iNum=1,nNum
drudeVB(:,1,iNum)=deg*drudeVB0(:,iNum)
drudeCB(:,1,iNum)=deg*drudeCB0(:,iNum)
enddo

dos_g=0._dp
drude_g=0._dp

vfermi=sqrt(3._dp)/2._dp*t0



do imu=1,nomega
fmu=real(imu-nshift)/real(nspacing)

do n=1,imu

!occCB(imu,1)=occCB(imu,1)+deg*dosCB0(n)/real(nspacing)*sqrt(3._dp)/2_dp*ndim/4_dp/Ac ! normalization to particle number
!occVB(imu,1)=occVB(imu,1)+deg*dosVB0(n)/real(nspacing)*sqrt(3._dp)/2_dp*ndim/4_dp/Ac ! normalization to particle number

occCB(imu,1)=occCB(imu,1)+deg*dosCB0(n)/real(nspacing)/Ac/A_n ! 10^12 cm^-2
occVB(imu,1)=occVB(imu,1)+deg*dosVB0(n)/real(nspacing)/Ac/A_n ! 10^12 cm^-2

enddo



do n=nshift+1,imu

ome=real(n-nshift)/real(nspacing)

occ_g(imu,1)=occ_g(imu,1)+degC*ome/vfermi**2/2._dp/pi/real(nspacing)/A_n

drude_g(imu,1)=drude_g(imu,1)+degC*ome/4._dp/pi

enddo


do iT=2,nT ! iT

T=ftemperature(iT)*T_in_ev

do n=nshift+1,nomega

ome=real(n-nshift)/real(nspacing)

occ_g(imu,iT)=occ_g(imu,iT)+degC*ome/vfermi**2/2._dp/pi/real(nspacing)*fermi_dist(ome-fmu,ftemperature(iT))/A_n

drude_g(imu,iT)=drude_g(imu,iT)+degC*exp((ome-fmu)/T)/(exp((ome-fmu)/T)+1._dp)**2/4._dp/pi*ome/T/real(nspacing)

enddo

enddo
do iNum=1,1
do iT=2,nT ! iT

T=ftemperature(iT)*T_in_ev ! Temperature in eV

!!$omp parallel do &
!!$omp private(n) &
!!$omp shared(imu,ome,fmu,T,dosCB0,dosVB0,drudeCB0,drudeVB0) &
!!$omp shared(nshift,nspacing,nstep) &
!!$omp reduction(+:drudeCB) &
!!$omp reduction(+:drudeVB) &
!!$omp reduction(+:occCB) &
!!$omp reduction(+:occVB)
do n=1,nomega,nstep

ome=real(n-nshift)/real(nspacing)

occCB(imu,iT)=occCB(imu,iT)+deg*dosCB0(n)/real(nspacing/nstep)*fermi_dist(ome-fmu,ftemperature(iT))/Ac/A_n
occVB(imu,iT)=occVB(imu,iT)+deg*dosVB0(n)/real(nspacing/nstep)*fermi_dist(ome-fmu,ftemperature(iT))/Ac/A_n

drudeCB(imu,iT,iNum)=drudeCB(imu,iT,iNum)+deg*exp((ome-fmu)/T)/(exp((ome-fmu)/T)+1._dp)**2*drudeCB0(n,iNum)/T/real(nspacing/nstep)
drudeVB(imu,iT,iNum)=drudeVB(imu,iT,iNum)+deg*exp((ome-fmu)/T)/(exp((ome-fmu)/T)+1._dp)**2*drudeVB0(n,iNum)/T/real(nspacing/nstep)

enddo
!!$omp end parallel do


enddo ! iT
enddo ! iNum

do iNum=2,nNum
do iT=2,nT ! iT

T=ftemperature(iT)*T_in_ev ! Temperature in eV

!!$omp parallel do &
!!$omp private(n) &
!!$omp shared(imu,ome,fmu,T,dosCB0,dosVB0,drudeCB0,drudeVB0) &
!!$omp shared(nshift,nspacing,nstep) &
!!$omp reduction(+:drudeCB) &
!!$omp reduction(+:drudeVB) &
!!$omp reduction(+:occCB) &
!!$omp reduction(+:occVB)
do n=1,nomega,nstep

ome=real(n-nshift)/real(nspacing)

drudeCB(imu,iT,iNum)=drudeCB(imu,iT,iNum)+deg*exp((ome-fmu)/T)/(exp((ome-fmu)/T)+1._dp)**2*drudeCB0(n,iNum)/T/real(nspacing/nstep)
drudeVB(imu,iT,iNum)=drudeVB(imu,iT,iNum)+deg*exp((ome-fmu)/T)/(exp((ome-fmu)/T)+1._dp)**2*drudeVB0(n,iNum)/T/real(nspacing/nstep)

enddo
!!$omp end parallel do


enddo ! iT
enddo ! iNum


enddo



! In units of e^2/hbar=1 and without spin and valley channel

THz=0.00413566553  !THz/eV

units=t0/THz*4._dp*pi**2

!write(*,*) units



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        do n=nshift,nomega

          write(101,fmt='(F16.8,3X)', advance="no") real(n-nshift)/real(nspacing)
          write(101,fmt='(F16.8,3X)', advance="no") deg*drudeVB0(n,1)/Ac
          write(101,fmt='(F16.8,3X)', advance="no") deg*drudeCB0(n,1)/Ac
          write(101,fmt='(F16.8,3X)') 2._dp/pi*real(n-nshift)/real(nspacing)

        end do

        close(101)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do iNum=1,nNum
     do iT=1,nT


        do n=nomega,2*nshift,-1

           write(110+iT+(iNum-1)*nNum,*) -occVB(n,iT),drudeVB(n,iT,iNum)/Ac,-occVB(n,iT)/A_n

    enddo

    do n=1,2*nshift-1

           write(110+iT+(iNum-1)*nNum,*) (occCB(n,iT)-occVB(2*nshift-n,iT)),(drudeCB(n,iT,iNum)+drudeVB(2*nshift-n,iT,iNum))/Ac,(occCB(n,iT)-occVB(2*nshift-n,iT))/A_n

        end do

    do n=2*nshift,nomega,1

           write(110+iT+(iNum-1)*nNum,*) occCB(n,iT),drudeCB(n,iT,iNum)/Ac,occCB(n,iT)/A_n

    enddo

        close(110+iT)

    enddo
    enddo

    do iT=1,nT

    if(iT.NE.1)then
    T=ftemperature(iT)*T_in_ev
    Lorentzian=T/pi*(0.0033_dp**2+T**2)
    else
    Lorentzian=1._dp
    endif

    enddo



    do iT=1,nT

        do n=nomega,2*nshift,-1

           write(120+iT,*) -occ_g(n,iT),drude_g(n,iT),-occ_g(n,iT)/A_n

    enddo

    do n=1,2*nshift-1

          write(120+iT,*) (occ_g(n,iT)-occ_g(2*nshift-n,iT)),(drude_g(n,iT)+drude_g(2*nshift-n,iT)),(occ_g(n,iT)-occ_g(2*nshift-n,iT))/A_n

        end do

    do n=2*nshift,nomega

           write(120+iT,*) occ_g(n,iT),drude_g(n,iT),occ_g(n,iT)/A_n

    enddo

        close(120+iT)

  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


iT=nT

        do n=1,nomega

           write(131,*) real(n-nshift)/real(nspacing),occVB(n,iT)/A_0

           write(132,*) real(n-nshift)/real(nspacing),occCB(n,iT)/A_0

        end do

    close(131)
    close(132)


end subroutine DrudeT


subroutine DrudeTemp(dos,drudeTot,drudeMag,drudeChi,nshift,nomega,nspacing,nNum)

 integer(dp), intent(in) :: nshift,nomega,nspacing,nNum
 real(dp)   , intent(in) :: dos(:),drudeTot(:),drudeMag(:),drudeChi(:)

 integer(dp) :: i,n,imu,iT,iNum
 real(dp) :: A_0,A_n,Ac,vfermi,T,fmu,ome,units,THz,Lorentzian

 integer(dp) , parameter :: nstep=1
 !integer(dp) , parameter :: nT=3
 !real(dp)    , parameter :: ftemperature(nT) =[0._dp,40._dp,170._dp]

 real(dp), parameter :: deg = 2._dp
 real(dp), parameter :: degC = 4._dp
 real(dp), parameter :: t0 = 2.7_dp
 real(dp), allocatable :: occ(:,:),occ_g(:,:)
 real(dp), allocatable :: drude(:,:,:),drude_g(:,:)



 allocate(occ(nomega,nT),occ_g(nomega,nT))
 allocate(drude(nomega,nT,nNum),drude_g(nomega,nT))

A_0=sqrt(3._dp)/2._dp*ndim/4._dp ! unit of density of one moire cell

A_n=2.46*2.46/10000._dp ! 10^-10/10000 = 10^-14  = 10^-12 (10^-2 m)= 10^-12 cm^-2

Ac=4._dp*pi*pi ! normalization

drude=0._dp

!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! Single-layer Graphene
!!!!!!!!!!!!!!!!!!!!!!

occ_g=0._dp
drude_g=0._dp


vfermi=sqrt(3._dp)/2._dp*t0

do imu=1,nomega

ome=real(imu-nshift)/real(nspacing)

occ_g(imu,1)=occ_g(imu,1)+degC*abs(ome)/vfermi**2/2._dp/pi/real(nspacing)/A_n

drude_g(imu,1)=drude_g(imu,1)+degC*abs(ome)/4._dp/pi

enddo


do iT=2,nT ! iT

T=ftemperature(iT)*T_in_ev


do imu=1,nomega
fmu=real(imu-nshift)/real(nspacing)
do n=1,nomega

ome=real(n-nshift)/real(nspacing)

occ_g(imu,iT)=occ_g(imu,iT)+degC*abs(ome)/vfermi**2/2._dp/pi/real(nspacing)*fermi_dist(ome-fmu,ftemperature(iT))/A_n

drude_g(imu,iT)=drude_g(imu,iT)+degC*exp((ome-fmu)/T)/(exp((ome-fmu)/T)+1._dp)**2/4._dp/pi*abs(ome)/T/real(nspacing)

enddo
enddo

enddo
!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Determine occupation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

occ=0._dp

do imu=2,nomega

!occ(imu,1)=occ(imu-1,1)+deg*dos(imu)/real(nspacing)*sqrt(3._dp)/2_dp*ndim/4_dp/Ac ! normalization to particle number

occ(imu,1)=occ(imu-1,1)+deg*dos(imu)/real(nspacing)/Ac/A_n ! 10^12 cm^-2

enddo

do imu=1,nomega
fmu=real(imu-nshift)/real(nspacing)

do iT=2,nT ! iT

T=ftemperature(iT)*T_in_ev ! Temperature in eV

!!$omp parallel do &
!!$omp private(n) &
!!$omp shared(imu,ome,fmu,T,dos) &
!!$omp shared(nshift,nspacing,nstep) &
!!$omp reduction(+:occ)
do n=1,nomega,nstep

ome=real(n-nshift)/real(nspacing)

occ(imu,iT)=occ(imu,iT)+deg*dos(n)/real(nspacing/nstep)*fermi_dist(ome-fmu,ftemperature(iT))/Ac/A_n

enddo
!!$omp end parallel do


enddo ! iT

enddo ! imu


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Determine Drude
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

drude(:,1,1)=deg*drudeTot(:)
drude(:,1,2)=deg*drudeMag(:)
drude(:,1,3)=deg*drudeChi(:)

do imu=1,nomega
fmu=real(imu-nshift)/real(nspacing)

do iNum=1,nNum
do iT=2,nT ! iT

T=ftemperature(iT)*T_in_ev ! Temperature in eV

!!$omp parallel do &
!!$omp private(n) &
!!$omp shared(imu,ome,fmu,T,drude,drude0) &
!!$omp shared(nshift,nspacing,nstep) &
!!$omp reduction(+:drude)
do n=1,nomega,nstep

ome=real(n-nshift)/real(nspacing)

drude(imu,iT,iNum)=drude(imu,iT,iNum)+1._dp/(exp(-(ome-fmu)/T)+1._dp)/(exp((ome-fmu)/T)+1._dp)*drude(n,1,iNum)/T/real(nspacing/nstep)

enddo
!!$omp end parallel do


enddo ! iT
enddo ! iNum


enddo ! imu



! In units of e^2/hbar=1 and without spin and valley channel

THz=0.00413566553  !THz/eV

units=t0/THz*4._dp*pi**2

!write(*,*) units



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do n=1,nomega

          write(101,fmt='(F16.8,3X)', advance="no") real(n-nshift)/real(nspacing)
          write(101,fmt='(F16.8,3X)', advance="no") drude(n,1,1)/Ac
          write(101,fmt='(F16.8,3X)') 2._dp/pi*abs(real(n-nshift))/real(nspacing)

        end do

        close(101)

        do n=1,nomega

           write(102,*) occ(n,1)-occ(nshift,1),deg*dos(n)*sqrt(3._dp)/2_dp*ndim/4_dp/Ac

    enddo

        close(102)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do iNum=1,nNum
     do iT=1,nT


        do n=1,nomega

           write(110+iT+10*(iNum-1),*) occ(n,iT)-occ(nshift,1),drude(n,iT,iNum)/Ac

    enddo

        close(110+iT+10*(iNum-1))

    enddo
    enddo

    do iT=1,nT

    if(iT.NE.1)then
    T=ftemperature(iT)*T_in_ev
    Lorentzian=T/pi*(0.0033_dp**2+T**2)
    else
    Lorentzian=1._dp
    endif

    enddo



    do iT=1,nT

        do n=1,nomega

           write(140+iT,*) occ_g(n,iT)-occ_g(n,1),drude_g(n,iT)

    enddo

        close(140+iT)

    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 do iT=1,nT

        do n=1,nomega

           write(150+iT,*) real(n-nshift)/real(nspacing),occ(n,iT)/A_0

        end do

    close(150+iT)

enddo

end subroutine DrudeTemp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine conductivityFullX(sigmaTot,sigmaMag,sigmaChi,vk,energia,pivotVx1,pivotVy1,pivotVx2,pivotVy2,pivotVxI,pivotVyI,nomega,nspacing,nshift,fmu,T,nbounds)

 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: fmu,T
 real(dp)   , intent(in) :: vk(:,:),energia(:,:)
 complex(dp), intent(in) :: pivotVx1(:,:,:),pivotVy1(:,:,:),pivotVx2(:,:,:),pivotVy2(:,:,:),pivotVxI(:,:,:),pivotVyI(:,:,:)
 real(dp)   , intent(out) :: sigmaTot(:),sigmaMag(:),sigmaChi(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE
 real(dp) :: energy(4),fermiDiff(4),weightTot(4),weightMag(4),weightChi(4)
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2),f1(3),f2(3),f3(3)
 real(dp) :: add,weight

  sigmaTot=0._dp
  sigmaMag=0._dp
  sigmaChi=0._dp


 do i_v = nbounds(1),nbounds(2)
   do i_c = nbounds(3),nbounds(4)

   !if(minval(fmu-energia(i_v,:)).GT.maxval(energia(i_c,:)-fmu))then
   !if(i_v.LT.i_c)then

   energy(:)=energia(i_c,:)-energia(i_v,:)

   !if(minval(energy).GT.0.00000001)then

   fermiDiff(:)=fermi_dist(energia(i_v,:)-fmu,T)-fermi_dist(energia(i_c,:)-fmu,T)

   weightTot(:)=real(conjg(pivotVx1(i_v,i_c,:)+pivotVx2(i_v,i_c,:)+pivotVxI(i_v,i_c,:))*(pivotVx1(i_v,i_c,:)+pivotVx2(i_v,i_c,:)+pivotVxI(i_v,i_c,:)))
   weightTot(:)=weightTot(:)/(energy(:)+0.0000000001_dp)*fermiDiff(:)

   weightMag(:)=real(conjg(pivotVx1(i_v,i_c,:)-pivotVx2(i_v,i_c,:))*(pivotVx1(i_v,i_c,:)-pivotVx2(i_v,i_c,:)))
   weightMag(:)=weightMag(:)/(energy(:)+0.0000000001_dp)*fermiDiff(:)

   weightChi(:)=real(conjg(pivotVx1(i_v,i_c,:))*(pivotVy2(i_v,i_c,:)+pivotVyI(i_v,i_c,:)))
   weightChi(:)=weightChi(:)/(energy(:)+0.0000000001_dp)*fermiDiff(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(:)=reshape([weightTot(1),weightMag(1),weightChi(1)],[3_dp])

   f2(:)=reshape([weightTot(2),weightMag(2),weightChi(2)],[3_dp])

   f3(:)=reshape([weightTot(3),weightMag(3),weightChi(3)],[3_dp])

   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*pi*real(iE-n3)/real(n2-n3)
         sigmaTot(iE)=sigmaTot(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         sigmaMag(iE)=sigmaMag(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
         sigmaChi(iE)=sigmaChi(iE)+add*(f3(3)*(n2-iE)/real(n2-n3)+f2(3)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*pi*real(n1-iE)/real(n1-n2)
        sigmaTot(iE)=sigmaTot(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        sigmaMag(iE)=sigmaMag(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
        sigmaChi(iE)=sigmaChi(iE)+add*(f2(3)*(n1-iE)/real(n1-n2)+f1(3)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(:)=reshape([weightTot(4),weightMag(4),weightChi(4)],[3_dp])

   f2(:)=reshape([weightTot(2),weightMag(2),weightChi(2)],[3_dp])

   f3(:)=reshape([weightTot(3),weightMag(3),weightChi(3)],[3_dp])

   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*pi*real(iE-n3)/real(n2-n3)
         sigmaTot(iE)=sigmaTot(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         sigmaMag(iE)=sigmaMag(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
         sigmaChi(iE)=sigmaChi(iE)+add*(f3(3)*(n2-iE)/real(n2-n3)+f2(3)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*pi*real(n1-iE)/real(n1-n2)
        sigmaTot(iE)=sigmaTot(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        sigmaMag(iE)=sigmaMag(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
        sigmaChi(iE)=sigmaChi(iE)+add*(f2(3)*(n1-iE)/real(n1-n2)+f1(3)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    !endif

    end do
    end do

end subroutine conductivityFullX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine conductivityFullY(sigmaTot,sigmaMag,sigmaChi,vk,energia,pivotVx1,pivotVy1,pivotVx2,pivotVy2,pivotVxI,pivotVyI,nomega,nspacing,nshift,fmu,T,nbounds)

 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: fmu,T
 real(dp)   , intent(in) :: vk(:,:),energia(:,:)
 complex(dp), intent(in) :: pivotVx1(:,:,:),pivotVy1(:,:,:),pivotVx2(:,:,:),pivotVy2(:,:,:),pivotVxI(:,:,:),pivotVyI(:,:,:)
 real(dp)   , intent(out) :: sigmaTot(:),sigmaMag(:),sigmaChi(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE
 real(dp) :: energy(4),fermiDiff(4),weightTot(4),weightMag(4),weightChi(4)
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2),f1(3),f2(3),f3(3)
 real(dp) :: add,weight

  sigmaTot=0._dp
  sigmaMag=0._dp
  sigmaChi=0._dp


 do i_v = nbounds(1),nbounds(2)
   do i_c = nbounds(3),nbounds(4)

   !if(minval(fmu-energia(i_v,:)).GT.maxval(energia(i_c,:)-fmu))then
   !if(i_v.LT.i_c)then

   energy(:)=energia(i_c,:)-energia(i_v,:)

   !if(minval(energy).GT.0.00000001)then

   fermiDiff(:)=fermi_dist(energia(i_v,:)-fmu,T)-fermi_dist(energia(i_c,:)-fmu,T)

   weightTot(:)=real(conjg(pivotVy1(i_v,i_c,:)+pivotVy2(i_v,i_c,:)+pivotVyI(i_v,i_c,:))*(pivotVy1(i_v,i_c,:)+pivotVy2(i_v,i_c,:)+pivotVyI(i_v,i_c,:)))
   weightTot(:)=weightTot(:)/(energy(:)+0.0000000001_dp)*fermiDiff(:)

   weightMag(:)=real(conjg(pivotVy1(i_v,i_c,:)-pivotVy2(i_v,i_c,:))*(pivotVy1(i_v,i_c,:)-pivotVy2(i_v,i_c,:)))
   weightMag(:)=weightMag(:)/(energy(:)+0.0000000001_dp)*fermiDiff(:)

   weightChi(:)=-real(conjg(pivotVy1(i_v,i_c,:))*(pivotVx2(i_v,i_c,:)+pivotVxI(i_v,i_c,:)))
   weightChi(:)=weightChi(:)/(energy(:)+0.0000000001_dp)*fermiDiff(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(:)=reshape([weightTot(1),weightMag(1),weightChi(1)],[3_dp])

   f2(:)=reshape([weightTot(2),weightMag(2),weightChi(2)],[3_dp])

   f3(:)=reshape([weightTot(3),weightMag(3),weightChi(3)],[3_dp])

   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*pi*real(iE-n3)/real(n2-n3)
         sigmaTot(iE)=sigmaTot(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         sigmaMag(iE)=sigmaMag(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
         sigmaChi(iE)=sigmaChi(iE)+add*(f3(3)*(n2-iE)/real(n2-n3)+f2(3)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*pi*real(n1-iE)/real(n1-n2)
        sigmaTot(iE)=sigmaTot(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        sigmaMag(iE)=sigmaMag(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
        sigmaChi(iE)=sigmaChi(iE)+add*(f2(3)*(n1-iE)/real(n1-n2)+f1(3)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(:)=reshape([weightTot(4),weightMag(4),weightChi(4)],[3_dp])

   f2(:)=reshape([weightTot(2),weightMag(2),weightChi(2)],[3_dp])

   f3(:)=reshape([weightTot(3),weightMag(3),weightChi(3)],[3_dp])

   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*pi*real(iE-n3)/real(n2-n3)
         sigmaTot(iE)=sigmaTot(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         sigmaMag(iE)=sigmaMag(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
         sigmaChi(iE)=sigmaChi(iE)+add*(f3(3)*(n2-iE)/real(n2-n3)+f2(3)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*pi*real(n1-iE)/real(n1-n2)
        sigmaTot(iE)=sigmaTot(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        sigmaMag(iE)=sigmaMag(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
        sigmaChi(iE)=sigmaChi(iE)+add*(f2(3)*(n1-iE)/real(n1-n2)+f1(3)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    !endif

    end do
    end do

end subroutine conductivityFullY

!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!

subroutine drudeFullX(Dtot,Dmag,Dxy,vk,energia,vx1,vy1,vx2,vy2,vxI,vyI,nomega,nspacing,nshift,efermi,nbounds)

 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: efermi,vk(:,:),energia(:,:),vx1(:,:),vy1(:,:),vx2(:,:),vy2(:,:),vxI(:,:),vyI(:,:)
 real(dp)   , intent(out) :: Dtot(:),Dmag(:),Dxy(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4)
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: f1(3),f2(3),f3(3),fvx1(4),fvx2(4),fvxI(4),fvy1(4),fvy2(4),fvyI(4),add,weight


   Dtot=0._dp
   Dmag=0._dp
   Dxy=0._dp

   do i_c = nbounds(1),nbounds(4)

   energy(:)=energia(i_c,:)-efermi

   fvx1(:)=vx1(i_c,:)

   fvx2(:)=vx2(i_c,:)

   fvxI(:)=vxI(i_c,:)

   fvy1(:)=vy1(i_c,:)

   fvy2(:)=vy2(i_c,:)

   fvyI(:)=vyI(i_c,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=(fvx1(1)+fvx2(1)+fvxI(1))*(fvx1(1)+fvx2(1)+fvxI(1))
   f2(1)=(fvx1(2)+fvx2(2)+fvxI(2))*(fvx1(2)+fvx2(2)+fvxI(2))
   f3(1)=(fvx1(3)+fvx2(3)+fvxI(3))*(fvx1(3)+fvx2(3)+fvxI(3))

   f1(2)=(fvx1(1)-fvx2(1))*(fvx1(1)-fvx2(1))
   f2(2)=(fvx1(2)-fvx2(2))*(fvx1(2)-fvx2(2))
   f3(2)=(fvx1(3)-fvx2(3))*(fvx1(3)-fvx2(3))

   f1(3)=fvx1(1)*(fvy2(1)+fvyI(1))
   f2(3)=fvx1(2)*(fvy2(2)+fvyI(2))
   f3(3)=fvx1(3)*(fvy2(3)+fvyI(3))




   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift


    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         Dtot(iE)=Dtot(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         Dmag(iE)=Dmag(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
         Dxy(iE)=Dxy(iE)+add*(f3(3)*(n2-iE)/real(n2-n3)+f2(3)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        Dtot(iE)=Dtot(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        Dmag(iE)=Dmag(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
        Dxy(iE)=Dxy(iE)+add*(f2(3)*(n1-iE)/real(n1-n2)+f1(3)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=(fvx1(4)+fvx2(4)+fvxI(4))*(fvx1(4)+fvx2(4)+fvxI(4))
   f2(1)=(fvx1(2)+fvx2(2)+fvxI(2))*(fvx1(2)+fvx2(2)+fvxI(2))
   f3(1)=(fvx1(3)+fvx2(3)+fvxI(3))*(fvx1(3)+fvx2(3)+fvxI(3))

   f1(2)=(fvx1(4)-fvx2(4))*(fvx1(4)-fvx2(4))
   f2(2)=(fvx1(2)-fvx2(2))*(fvx1(2)-fvx2(2))
   f3(2)=(fvx1(3)-fvx2(3))*(fvx1(3)-fvx2(3))

   f1(3)=fvx1(4)*(fvy2(4)+fvyI(4))
   f2(3)=fvx1(2)*(fvy2(2)+fvyI(2))
   f3(3)=fvx1(3)*(fvy2(3)+fvyI(3))

   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         Dtot(iE)=Dtot(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         Dmag(iE)=Dmag(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
         Dxy(iE)=Dxy(iE)+add*(f3(3)*(n2-iE)/real(n2-n3)+f2(3)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        Dtot(iE)=Dtot(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        Dmag(iE)=Dmag(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
        Dxy(iE)=Dxy(iE)+add*(f2(3)*(n1-iE)/real(n1-n2)+f1(3)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    end do

    !write(*,*) 'TestDos',maxval(dos)

end subroutine drudeFullX

!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!

subroutine drudeFullY(Dtot,Dmag,Dxy,vk,energia,vx1,vy1,vx2,vy2,vxI,vyI,nomega,nspacing,nshift,efermi,nbounds)

 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: efermi,vk(:,:),energia(:,:),vx1(:,:),vy1(:,:),vx2(:,:),vy2(:,:),vxI(:,:),vyI(:,:)
 real(dp)   , intent(out) :: Dtot(:),Dmag(:),Dxy(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4)
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: f1(3),f2(3),f3(3),fvx1(4),fvx2(4),fvxI(4),fvy1(4),fvy2(4),fvyI(4),add,weight


   Dtot=0._dp
   Dmag=0._dp
   Dxy=0._dp

   do i_c = nbounds(1),nbounds(4)

   energy(:)=energia(i_c,:)-efermi

   fvx1(:)=vx1(i_c,:)

   fvx2(:)=vx2(i_c,:)

   fvxI(:)=vxI(i_c,:)

   fvy1(:)=vy1(i_c,:)

   fvy2(:)=vy2(i_c,:)

   fvyI(:)=vyI(i_c,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=(fvy1(1)+fvy2(1)+fvyI(1))*(fvy1(1)+fvy2(1)+fvyI(1))
   f2(1)=(fvy1(2)+fvy2(2)+fvyI(2))*(fvy1(2)+fvy2(2)+fvyI(2))
   f3(1)=(fvy1(3)+fvy2(3)+fvyI(3))*(fvy1(3)+fvy2(3)+fvyI(3))

   f1(2)=(fvy1(1)-fvy2(1))*(fvy1(1)-fvy2(1))
   f2(2)=(fvy1(2)-fvy2(2))*(fvy1(2)-fvy2(2))
   f3(2)=(fvy1(3)-fvy2(3))*(fvy1(3)-fvy2(3))

   f1(3)=-(fvx2(1)+fvxI(1))*fvy1(1)
   f2(3)=-(fvx2(2)+fvxI(2))*fvy1(2)
   f3(3)=-(fvx2(3)+fvxI(3))*fvy1(3)




   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift


    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         Dtot(iE)=Dtot(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         Dmag(iE)=Dmag(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
         Dxy(iE)=Dxy(iE)+add*(f3(3)*(n2-iE)/real(n2-n3)+f2(3)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        Dtot(iE)=Dtot(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        Dmag(iE)=Dmag(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
        Dxy(iE)=Dxy(iE)+add*(f2(3)*(n1-iE)/real(n1-n2)+f1(3)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=(fvy1(4)+fvy2(4)+fvyI(4))*(fvy1(4)+fvy2(4)+fvyI(4))
   f2(1)=(fvy1(2)+fvy2(2)+fvyI(2))*(fvy1(2)+fvy2(2)+fvyI(2))
   f3(1)=(fvy1(3)+fvy2(3)+fvyI(3))*(fvy1(3)+fvy2(3)+fvyI(3))

   f1(2)=(fvy1(4)-fvy2(4))*(fvy1(4)-fvy2(4))
   f2(2)=(fvy1(2)-fvy2(2))*(fvy1(2)-fvy2(2))
   f3(2)=(fvy1(3)-fvy2(3))*(fvy1(3)-fvy2(3))

   f1(3)=-(fvx2(4)+fvxI(4))*fvy1(4)
   f2(3)=-(fvx2(2)+fvxI(2))*fvy1(2)
   f3(3)=-(fvx2(3)+fvxI(3))*fvy1(3)

   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         Dtot(iE)=Dtot(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
         Dmag(iE)=Dmag(iE)+add*(f3(2)*(n2-iE)/real(n2-n3)+f2(2)*(iE-n3)/real(n2-n3))
         Dxy(iE)=Dxy(iE)+add*(f3(3)*(n2-iE)/real(n2-n3)+f2(3)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        Dtot(iE)=Dtot(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
        Dmag(iE)=Dmag(iE)+add*(f2(2)*(n1-iE)/real(n1-n2)+f1(2)*(iE-n2)/real(n1-n2))
        Dxy(iE)=Dxy(iE)+add*(f2(3)*(n1-iE)/real(n1-n2)+f1(3)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    end do

    !write(*,*) 'TestDos',maxval(dos)

end subroutine drudeFullY

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%^M
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%^M

subroutine dos_F(dos,vk,energia,nomega,nspacing,nshift,fmu,nbounds)


 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: vk(:,:),energia(:,:),fmu
 real(dp)   , intent(out) :: dos(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4)
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: add,weight


   dos=0._dp

   do i_c = nbounds(1),nbounds(4)

   energy(:)=energia(i_c,:)-fmu


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   call sort(energy1,energy2,energy3,vk1,vk2,vk3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift


    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         dos(iE)=dos(iE)+add
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        dos(iE)=dos(iE)+add
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   call sort(energy1,energy2,energy3,vk1,vk2,vk3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         dos(iE)=dos(iE)+add
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        dos(iE)=dos(iE)+add
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    end do

    !write(*,*) 'TestDos',maxval(dos)

end subroutine dos_F

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine KramersKronig(ReSigma,ImSigma,nKK,nomega,nspacing,nshift,nstep,nstepI)

 integer(dp), intent(in) :: nKK,nomega,nspacing,nshift,nstep,nstepI
 real(dp)   , intent(in) :: ImSigma(:)
 real(dp)   , intent(out) :: ReSigma(:)

 integer(dp) :: i,k
 integer(dp) :: nomegaMax
 real(dp) :: omega(nomega)

 ReSigma=0._dp
do i=nshift,nomega
omega(i)=real(i-nshift,dp)/real(nspacing,dp)
enddo

nomegaMax=nomega
do while(nomegaMax >= nshift)
    if(ImSigma(nomegaMax) >= 0.0000000001_dp) exit
    nomegaMax = nomegaMax-1
enddo
if(nomegaMax < nshift) return

write(*,*) 'nomegaMax',nomegaMax,nomega

    do k=nshift,nKK,nstep
!$omp parallel do &
!$omp private(i) &
!$omp shared(k,ImSigma,omega,nstepI) &
!$omp reduction(+:ReSigma)
    do i=nshift,k-nstepI,nstepI
      ReSigma(k)=ReSigma(k)+omega(i)*omega(i)*ImSigma(i)/(omega(k)**2-omega(i)**2)
    enddo
!$omp end parallel do

!$omp parallel do &
!$omp private(i) &
!$omp shared(k,ImSigma,omega,nstepI) &
!$omp reduction(+:ReSigma)
    do i=k+nstepI,nomegaMax,nstepI
      ReSigma(k)=ReSigma(k)+omega(i)*omega(i)*ImSigma(i)/(omega(k)**2-omega(i)**2)
    enddo
!$omp end parallel do

    enddo

    ReSigma=ReSigma/real(nspacing/nstepI,dp)

end subroutine KramersKronig

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine conductivityS(sigmaAdd,vk,energia,pivotVx,pivotVy,nomega,nspacing,nshift,fmu,T,nbounds)

 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: fmu,T
 real(dp)   , intent(in) :: vk(:,:),energia(:,:)
 complex(dp), intent(in) :: pivotVx(:,:,:),pivotVy(:,:,:)
 real(dp)   , intent(out) :: sigmaAdd(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE
 real(dp) :: energy(4),fermiDiff(4),weightAdd(4)
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2),f1(3),f2(3),f3(3)
 real(dp) :: add,weight

 sigmaAdd=0._dp


 do i_v = nbounds(1),nbounds(2)
   do i_c = nbounds(3),nbounds(4)

   !if(minval(fmu-energia(i_v,:)).GT.maxval(energia(i_c,:)-fmu))then
   !if(i_v.LT.i_c)then
   !if(i_v.NE.i_c)then

   energy(:)=energia(i_c,:)-energia(i_v,:)

   !if(minval(energy).GT.0.00000001)then

   fermiDiff(:)=fermi_dist(energia(i_v,:)-fmu,T)-fermi_dist(energia(i_c,:)-fmu,T)

   weightAdd(:)=real(conjg(pivotVx(i_v,i_c,:))*pivotVx(i_v,i_c,:)+conjg(pivotVy(i_v,i_c,:)) &
   *pivotVy(i_v,i_c,:))/(energy(:)+0.0000000001_dp)*fermiDiff(:)/2._dp

   if(minval(weightAdd).LT.0.)then
           write(*,*) 'Alarm',weightAdd(:)
           write(*,*) i_v,i_c,energy
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate first triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

!   f1(:)=reshape([weightAdd(1),weightMinus(1),weightX1Y2(1)],[3_dp])

!   f2(:)=reshape([weightAdd(2),weightMinus(2),weightX1Y2(2)],[3_dp])

!   f3(:)=reshape([weightAdd(3),weightMinus(3),weightX1Y2(3)],[3_dp])

   f1(1)=weightAdd(1)

   f2(1)=weightAdd(2)

   f3(1)=weightAdd(3)

   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then

    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*pi*real(iE-n3)/real(n2-n3)
         sigmaAdd(iE)=sigmaAdd(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*pi*real(n1-iE)/real(n1-n2)
        sigmaAdd(iE)=sigmaAdd(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  Calculate second triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)

   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)

   f1(1)=weightAdd(4)

   f2(1)=weightAdd(2)

   f3(1)=weightAdd(3)

   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

       weight=abs(((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))/(energy1-energy3))

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift

    if(n1.LE.nomega) then
        if(n3.GE.1) then


    if(n3.LT.n2)then

    do iE=n3+1,n2
         add=weight*pi*real(iE-n3)/real(n2-n3)
         sigmaAdd(iE)=sigmaAdd(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
    end do

    endif

    if(n2.LT.n1)then

    do iE=n2+1,n1-1
        add=weight*pi*real(n1-iE)/real(n1-n2)
        sigmaAdd(iE)=sigmaAdd(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
    end do

    endif

    endif ! n3.GE.1
    endif ! n1.LE.nomega

    !endif

    end do
    end do

end subroutine conductivityS



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine plotBandsK(bandsTB,fmu,npointsBack,ndimEV,numk,ncount,vq1,vq12)

 integer(dp), intent(in)    :: ndimEV,numk,ncount
 integer(dp), intent(in)    :: npointsBack(numk,numk)
 real(dp)   , intent(in)    :: fmu,vq1(2), vq12(2), bandsTB(ndimEV,numk+1,numk+1)

 integer(dp) :: i,ivk,ivk1,ivk2
 real(dp) :: vK1x,vK1y,vK2x,vK2y,vM1x,vM1y,vM2x,vM2y
 real(dp) :: a11,a12,a21,a22,b11,b12,b21,b22,vkx,vky,det
 real(dp) :: aGM, aKG, aMK, aT

 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 vK1x=(vq1(1)+vq12(1))/3.0_dp
 vK1y=(vq1(2)+vq12(2))/3.0_dp

 vK2x=2.0_dp*(vq1(1)+vq12(1))/3.0_dp
 vK2y=2.0_dp*(vq1(2)+vq12(2))/3.0_dp

 vM1x=vq1(1)/2.0_dp
 vM1y=vq1(2)/2.0_dp

 vM2x=vq12(1)/2.0_dp
 vM2y=vq12(2)/2.0_dp


 aKG=sqrt(vK1x**2+vK1y**2)
 aGM=sqrt(vM1x**2+vM1y**2)
 aMK=sqrt((vK1x-vM1x)**2+(vK1y-vM1y)**2)
 aT=aKG+aGM+aMK

 a11=vq1(1)/real(numk,dp)
 a12=vq12(1)/real(numk,dp)
 a21=vq1(2)/real(numk,dp)
 a22=vq12(2)/real(numk,dp)

 det=a11*a22-a12*a21
 b11=a22/det
 b12=-a12/det
 b21=-a21/det
 b22=a11/det

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

 do ivk2=numk/3,0,-1
     ivk1=ivk2
     !write(*,*) (numk/3-ivk2)*aKG/real(numk/3,dp)/aT,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
     write(99,fmt='(F16.8,3X)', advance="no") (numk/3-ivk2)*aKG/real(numk/3,dp)/aT
     do i=1,ndimEV-1
     write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,ivk1+1,ivk2+1)-fmu
     enddo
     write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,ivk1+1,ivk2+1)-fmu
 enddo
 do ivk2=1,numk/2
     ivk1=0
     !write(*,*) (aKG+(ivk2)*aGM/real(numk/2,dp))/aT,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
     write(99,fmt='(F16.8,3X)', advance="no") (aKG+(ivk2)*aGM/real(numk/2,dp))/aT
     do i=1,ndimEV-1
     write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,ivk1+1,ivk2+1)-fmu
     enddo
     write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,ivk1+1,ivk2+1)-fmu
 enddo
 do ivk=1,numk/6
     vkx=vM1x+(vK1x-vM1x)*ivk/real(numk/6,dp)
     vky=vM1y+(vK1y-vM1y)*ivk/real(numk/6,dp)
     ivk2=nint(b11*vkx+b12*vky)
     ivk1=nint(b21*vkx+b22*vky)
     !write(*,*) (aKG+aGM+(ivk)*aMK/real(numk/6,dp))/aT,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
     write(99,fmt='(F16.8,3X)', advance="no") (aKG+aGM+(ivk)*aMK/real(numk/6,dp))/aT
     do i=1,ndimEV-1
     write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,ivk1+1,ivk2+1)-fmu
     enddo
     write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,ivk1+1,ivk2+1)-fmu
 enddo

 close(99)

end subroutine plotBandsK



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine plotBandsDensity_dB(bandsTB,npointsBZ,ndimEV,numk,nband,ncount,nNP)

 integer(dp), intent(in)    :: ndimEV,numk,nband,ncount,nNP
 integer(dp), intent(in)    :: npointsBZ(ncount,8)
 real(dp), intent(in)    :: bandsTB(ndimEV,ncount)

 integer(dp) :: i,icount,ivk1,ivk2

 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 do icount=1,ncount

     ivk1=npointsBZ(icount,1)
     ivk2=npointsBZ(icount,2)

     write(98,fmt='(F16.8,3X)', advance="no") real(ivk1,dp)
     write(98,fmt='(F16.8,3X)', advance="no") real(ivk2,dp)
     write(98,fmt='(F16.8,3X)') bandsTB(nNP+nband,icount)

enddo

end subroutine plotBandsDensity_dB


end module Hamiltonian
