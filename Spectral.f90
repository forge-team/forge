! Compatibility routines grouped by the same responsibilities as FORGE dev.
module Spectral
use ResponseSetup
use lapack_routines
use ResponseGeometry
use ResponseTightBinding
use ResponseHartreeFock
implicit none
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



! BEGIN imported local-DOS kernels
subroutine LocalDOS_CB(ndim,ldos,vk,energia,vPsi,nomega,nspacing,nshift,efermi,nbounds)
 integer(dp), intent(in) :: ndim,nomega,nspacing,nshift,nbounds(:)
 real(dp)   , intent(in) :: efermi
 real(dp)   , intent(in) :: vk(:,:),energia(:,:),vPsi(:,:)
 real(dp)   , intent(out) :: ldos(:)
 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4)
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: f1(1),f2(1),f3(1)
 real(dp) :: add,weight
   ldos=0._dp
   do i_c = nbounds(3),nbounds(4)
   energy(:)=energia(i_c,:)-efermi
   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)
   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)
   f1(1)=vPsi(i_c,1)
   f2(1)=vPsi(i_c,2)
   f3(1)=vPsi(i_c,3)
   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift
    ! An unresolved/flat triangle carries its integrated weight in one bin.
    weight=abs((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))
    if(n1 == n3) then
        if(n1 >= 1 .and. n1 <= nomega) ldos(n1)=ldos(n1)+weight*nspacing*(f1(1)+f2(1)+f3(1))/6.0_dp
    else
        weight=weight/abs(energy1-energy3)
    endif
    if(n1.LE.nomega) then
        if(n3.GE.1) then
    if(n3.LT.n2)then
    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         ldos(iE)=ldos(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
    end do
    endif
    if(n2.LT.n1)then
    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        ldos(iE)=ldos(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
    end do
    endif
    endif ! n3.GE.1
    endif ! n1.LE.nomega
   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)
   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)
   f1(1)=vPsi(i_c,4)
   f2(1)=vPsi(i_c,2)
   f3(1)=vPsi(i_c,3)
   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift
    ! An unresolved/flat triangle carries its integrated weight in one bin.
    weight=abs((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))
    if(n1 == n3) then
        if(n1 >= 1 .and. n1 <= nomega) ldos(n1)=ldos(n1)+weight*nspacing*(f1(1)+f2(1)+f3(1))/6.0_dp
    else
        weight=weight/abs(energy1-energy3)
    endif
    if(n1.LE.nomega) then
        if(n3.GE.1) then
    if(n3.LT.n2)then
    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         ldos(iE)=ldos(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
    end do
    endif
    if(n2.LT.n1)then
    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        ldos(iE)=ldos(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
    end do
    endif
    endif ! n3.GE.1
    endif ! n1.LE.nomega
    enddo
end subroutine LocalDOS_CB

subroutine LocalDOS_VB(ndim,ldos,vk,energia,vPsi,nomega,nspacing,nshift,efermi,nbounds)
 integer(dp), intent(in) :: ndim,nomega,nspacing,nshift,nbounds(:)
 real(dp)   , intent(in) :: efermi
 real(dp)   , intent(in) :: vk(:,:),energia(:,:),vPsi(:,:)
 real(dp)   , intent(out) :: ldos(:)
 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4)
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: f1(1),f2(1),f3(1)
 real(dp) :: add,weight
   ldos=0._dp
   do i_v = nbounds(1),nbounds(2)
   energy(:)=efermi-energia(i_v,:)
   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)
   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)
   f1(1)=vPsi(i_v,1)
   f2(1)=vPsi(i_v,2)
   f3(1)=vPsi(i_v,3)
   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift
    ! An unresolved/flat triangle carries its integrated weight in one bin.
    weight=abs((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))
    if(n1 == n3) then
        if(n1 >= 1 .and. n1 <= nomega) ldos(n1)=ldos(n1)+weight*nspacing*(f1(1)+f2(1)+f3(1))/6.0_dp
    else
        weight=weight/abs(energy1-energy3)
    endif
    if(n1.LE.nomega) then
        if(n3.GE.1) then
    if(n3.LT.n2)then
    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         ldos(iE)=ldos(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
    end do
    endif
    if(n2.LT.n1)then
    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        ldos(iE)=ldos(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
    end do
    endif
    endif ! n3.GE.1
    endif ! n1.LE.nomega
   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)
   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)
   f1(1)=vPsi(i_v,4)
   f2(1)=vPsi(i_v,2)
   f3(1)=vPsi(i_v,3)
   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift
    ! An unresolved/flat triangle carries its integrated weight in one bin.
    weight=abs((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))
    if(n1 == n3) then
        if(n1 >= 1 .and. n1 <= nomega) ldos(n1)=ldos(n1)+weight*nspacing*(f1(1)+f2(1)+f3(1))/6.0_dp
    else
        weight=weight/abs(energy1-energy3)
    endif
    if(n1.LE.nomega) then
        if(n3.GE.1) then
    if(n3.LT.n2)then
    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         ldos(iE)=ldos(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
    end do
    endif
    if(n2.LT.n1)then
    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        ldos(iE)=ldos(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
    end do
    endif
    endif ! n3.GE.1
    endif ! n1.LE.nomega
    end do
end subroutine LocalDOS_VB
! END imported local-DOS kernels

end module Spectral
