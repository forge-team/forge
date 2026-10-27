module TightBinding

use Setup
use lapack_routines

implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine TightBindingHamiltonian(zH,Delta,Coords,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,vk,tn,NearestNeighborsUC,NearestNeighborsT)
 
 integer(dp), intent(in)    :: ndim,numNeighborCells
 integer(dp), intent(in)    :: nUnitCell_1(numNeighborCells),nUnitCell_2(numNeighborCells),NearestNeighborsUC(ndim,3),NearestNeighborsT(ndim,3)
 real(dp)   , intent(in)    :: tn(2,3), vk(2), Coords(ndim,3), Delta
 complex(dp), intent(out) :: zH(ndim,ndim)
 
 integer(dp) :: i, j, n1,n2,icount
 real(dp)    :: rij, Coords_diff, tm(0:6,2)
 complex(dp) :: zi_vk_tn,zphase

  tm(0,:) = [0.0_dp, 0.0_dp]
  tm(1,:) = tn(:,1)
  tm(2,:) = tn(:,2)
  tm(3,:) = tn(:,3)
  tm(4,:) = -tn(:,1)
  tm(5,:) = -tn(:,2)
  tm(6,:) = -tn(:,3)
 
  zH = cmplx(0.0_dp,0.0_dp,dp)

  do i = 1 , ndim
    do j = 1 , i
   
      zphase=exp(cmplx(0.0_dp,-dot_product(vk,Coords(i,1:2)-Coords(j,1:2))*real(fphase,dp),dp))
     
      zH(i,j) = zH(i,j) + fTB(TBFunction,Coords(i,:),Coords(j,:),i,j,NearestNeighborsUC,NearestNeighborsT,Coords,tm)*zphase
      do icount=2,numNeighborCells
        n1=nUnitCell_1(icount)
        n2=nUnitCell_2(icount)
            
        zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
        zH(i,j)=zH(i,j)+fTB(TBFunction,Coords(i,:),[Coords(j,1:2)-n1*tn(:,1)-n2*tn(:,2),Coords(j,3)],i,j,NearestNeighborsUC,NearestNeighborsT,Coords,tm)*exp(-zi_vk_tn)*zphase
                
        zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
        zH(i,j)=zH(i,j)+fTB(TBFunction,Coords(i,:),[Coords(j,1:2)+n1*tn(:,1)+n2*tn(:,2),Coords(j,3)],i,j,NearestNeighborsUC,NearestNeighborsT,Coords,tm)*exp(-zi_vk_tn)*zphase

      enddo

    end do
  end do
  if(nlayers.GT.1)then
  do n1 = 0,nlayers-1
    do i=1,ndim/nlayers
      zH(i+ n1*ndim/nlayers,i+ n1*ndim/nlayers) = zH(i+ n1*ndim/nlayers,i+ n1*ndim/nlayers) + (n1 - real(nlayers-1)/2)*Delta/real(nlayers-1)
    enddo
  enddo
  endif
 
end subroutine TightBindingHamiltonian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine ValleyPhase(zV,Coords,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,vk,tn,RotMatrix)
 
  integer(dp), intent(in)    :: ndim,numNeighborCells
  integer(dp), intent(in)    :: nUnitCell_1(numNeighborCells),nUnitCell_2(numNeighborCells)
  real(dp)   , intent(in)    :: tn(2,3), vk(2), Coords(ndim,3), RotMatrix(2,2)
  complex(dp), intent(out) :: zV(ndim,ndim)
  
  integer(dp) :: i, j, n1,n2,icount,nlayer
  real(dp)    :: ri(3), rj(3)
  complex(dp) :: zi_vk_tn,zphase
 
  zV = cmplx(0.0_dp,0.0_dp,dp)
  
  do nlayer=1,nlayers

    do i = ndim/nlayers*(nlayer-1) + 1 ,ndim/nlayers*(nlayer-1) + ndim/nlayers/2
      do j = ndim/nlayers*(nlayer-1) + 1 ,ndim/nlayers*(nlayer-1) + ndim/nlayers/2

        ri =  Coords(i,:)
        rj =  Coords(j,:)
      
        zphase=exp(cmplx(0.0_dp,-dot_product(vk,Coords(i,1:2)-Coords(j,1:2))*real(fphase,dp),dp))
        
        if (RotateLayers(nlayer).eq.-1)then
          zV(i,j) = zV(i,j) + fValley(ri,rj,transpose(RotMatrix))*zphase
        else if (RotateLayers(nlayer).eq.1)then
          zV(i,j) = zV(i,j) + fValley(ri,rj,RotMatrix)*zphase
        endif

        do icount=2,numNeighborCells
          n1=nUnitCell_1(icount)
          n2=nUnitCell_2(icount)
        
          ri =  Coords(i,:)
          rj =  [Coords(j,1:2) - n1*tn(:,1) - n2*tn(:,2), Coords(j,3)]
          
          zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
          
          if (RotateLayers(nlayer).eq.-1)then
            zV(i,j) = zV(i,j) + fValley(ri,rj,transpose(RotMatrix))*zphase*exp(-zi_vk_tn)
          else if (RotateLayers(nlayer).eq.1)then
            zV(i,j) = zV(i,j) + fValley(ri,rj,RotMatrix)*zphase*exp(-zi_vk_tn)
          endif
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ri =  [Coords(i,1:2) - n1*tn(:,1) - n2*tn(:,2), Coords(i,3)]
          rj =  Coords(j,:) 
          
          zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
          
          if (RotateLayers(nlayer).eq.-1)then
            zV(i,j) = zV(i,j) + fValley(ri,rj,transpose(RotMatrix))*zphase*exp(-zi_vk_tn)
          else if (RotateLayers(nlayer).eq.1)then
            zV(i,j) = zV(i,j) + fValley(ri,rj,RotMatrix)*zphase*exp(-zi_vk_tn)
          endif           
        enddo
  
      end do
    end do
  
    do i = ndim/nlayers*(nlayer-1) + ndim/nlayers/2 + 1 ,ndim/nlayers*(nlayer-1) + ndim/nlayers
      do j = ndim/nlayers*(nlayer-1) + ndim/nlayers/2 + 1 ,ndim/nlayers*(nlayer-1) + ndim/nlayers

        ri =  Coords(i,:)
        rj =  Coords(j,:)
      
        zphase=exp(cmplx(0.0_dp,-dot_product(vk,Coords(i,1:2)-Coords(j,1:2))*real(fphase,dp),dp))
        
        if (RotateLayers(nlayer).eq.-1)then
          zV(i,j) = zV(i,j) + fValley(ri,rj,transpose(RotMatrix))*zphase
        else if (RotateLayers(nlayer).eq.1)then
          zV(i,j) = zV(i,j) + fValley(ri,rj,RotMatrix)*zphase
        endif

        do icount=2,numNeighborCells
          n1=nUnitCell_1(icount)
          n2=nUnitCell_2(icount)
        
          ri =  Coords(i,:)
          rj =  [Coords(j,1:2) - n1*tn(:,1) - n2*tn(:,2), Coords(j,3)]
          
          zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
          
          if (RotateLayers(nlayer).eq.-1)then
            zV(i,j) = zV(i,j) + fValley(ri,rj,transpose(RotMatrix))*zphase*exp(-zi_vk_tn)
          else if (RotateLayers(nlayer).eq.1)then
            zV(i,j) = zV(i,j) + fValley(ri,rj,RotMatrix)*zphase*exp(-zi_vk_tn)
          endif         
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ri =  [Coords(i,1:2) - n1*tn(:,1) - n2*tn(:,2), Coords(i,3)]
          rj =  Coords(j,:) 
          
          zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
          
          if (RotateLayers(nlayer).eq.-1)then
            zV(i,j) = zV(i,j) + fValley(ri,rj,transpose(RotMatrix))*zphase*exp(-zi_vk_tn)
          else if (RotateLayers(nlayer).eq.1)then
            zV(i,j) = zV(i,j) + fValley(ri,rj,RotMatrix)*zphase*exp(-zi_vk_tn)
          endif           
        enddo
  
      end do
    end do

  enddo
   
 end subroutine ValleyPhase

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine ValleyTransform(zEigenvectors,Coords,MomentaValues,nMomentaComponents,&
  numb,ndim,numk,Nk,numNeighborCells,nUnitCell_1,nUnitCell_2,tn,RotMatrix)

  integer(dp), intent(in)    :: ndim, numk, Nk, numb, numNeighborCells
  integer(dp), intent(in)    :: nMomentaComponents(Nk,2), nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells)
  real(dp)   , intent(in)    :: MomentaValues(numk,numk,2), Coords(ndim,3), tn(2,3),RotMatrix(2,2)
  complex(dp), intent(inout) :: zEigenvectors(ndim,numb,Nk)

  integer(dp) :: ivk1, ivk2, icount
  real(dp)    :: vk(2)
  complex(dp) :: zvalley(ndim,ndim), zvalleyproj(numb,numb)
  complex(dp) :: ztemp(ndim,numb)

  !$omp parallel do &
  !$omp private(icount,vk,ivk1,ivk2,zvalley,zvalleyproj,ztemp) &
  !$omp shared(Nk,nMomentaComponents,MomentaValues,Coords,ndim,numNeighborCells,nUnitCell_2,nUnitCell_1,tn,RotMatrix,numb,zEigenvectors)
  do icount=1,Nk
      ivk1=nMomentaComponents(icount,1)
      ivk2=nMomentaComponents(icount,2)
      vk(:)=MomentaValues(ivk1+1,ivk2+1,:)
      
      write(*,*) icount, Nk

      !write(*,*) 'Create valley matrix at vk'
      zvalley(:,:) = cmplx(0.0_dp,0.0_dp,dp)
      call ValleyPhase(zvalley,Coords,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,vk,tn,RotMatrix)

      zvalleyproj(:,:) = cmplx(0.0_dp,0.0_dp,dp)
      !write(*,*) 'Project valley matrix'
      call matrix_elements(zvalley,zEigenvectors(:,:,icount),zvalleyproj)
      
      !write(*,*) 'Unitarize valley matrix'
      call unitarize(zvalleyproj)
      
      !write(*,*) 'Compute tau_z*Eigenvectors'
      call matrix_mul(zEigenvectors(:,:,icount),zvalleyproj,ztemp(:,:))
      zEigenvectors(:,:,icount) = ztemp(:,:)
      
  enddo
  !$omp end parallel do

end subroutine ValleyTransform


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fTB(TBFunction,r_i,r_j,i,j,NearestNeighborsUC,NearestNeighborsT,Coords,tn)

  integer(dp), intent(in) :: TBFunction                   ! 1/2: Slater-Koster/Kaxiras hopping (TBFunction)
  integer(dp), intent(in) :: i,j, NearestNeighborsUC(ndim,3), NearestNeighborsT(ndim,3)
  real(dp), intent(in) :: r_i(3),r_j(3),tn(0:6,2), Coords(ndim,3)

  ! TBFunction=1: Slater-Koster
  real(dp), parameter :: r0 = 0.184_dp                   ! tight-binding decay constant
  real(dp), parameter :: d0 = 1.35772_dp        ! reference interlayer distance in units of a

  ! TBFunction=2: Wannier, intralayer hopping
  real(dp), parameter :: tparallel0 = -18.4295_dp
  real(dp), parameter :: al0 = 1.2771_dp
  real(dp), parameter :: bet0 = 2.3934_dp
  real(dp), parameter :: tparallel1 = -3.7183_dp
  real(dp), parameter :: al1 = 6.2194_dp
  real(dp), parameter :: y1 = 0.9071_dp
  ! TBFunction=2: Wannier, interlayer hopping
  real(dp), parameter :: lambda0 = 0.3155_dp
  real(dp), parameter :: lambda3 = -1*0.0688_dp
  real(dp), parameter :: lambda6 = -1*0.0083_dp
  real(dp), parameter :: xi0 = 1.7543_dp
  real(dp), parameter :: xi3 = 3.4692_dp
  real(dp), parameter :: xi6 = 2.8764_dp
  real(dp), parameter :: y3 = 0.5212_dp
  real(dp), parameter :: y6 = 1.5206_dp
  real(dp), parameter :: kappa0 = 2.0010_dp
  real(dp), parameter :: kappa6 = 1.5731_dp

  integer(dp) :: nni(3), nti(3), nnj(3), ntj(3)
  real(dp) :: theta12_1, theta12_2, theta12_3
  real(dp) :: theta21_1, theta21_2, theta21_3
  real(dp) :: V0,V3,V6,r,cs

  real(dp) :: fTB

  fTB = 0.0_dp

  select case(TBFunction)

  case(1)

    r = norm2(r_i-r_j)
    cs = (r_i(3)-r_j(3))/norm2(r_i-r_j)

    if(r.GT.0.001_dp)then
      fTB = -2.7_dp*exp((a0-r)/r0)*(1-cs**2) + 0.48_dp*exp((d0-r)/r0)*cs**2
    endif

  case(2)

    if(abs(r_i(3) - r_j(3)).gt.(tz*0.1))then

      nni = NearestNeighborsUC(i,:)
      nnj = NearestNeighborsUC(j,:)
      nti = NearestNeighborsT(i,:)
      ntj = NearestNeighborsT(j,:)

      theta21_1 = atan2(r_i(2)-r_j(2),r_i(1)-r_j(1))&
                  - atan2(Coords(nnj(1),2) - tn(ntj(1),2)-Coords(j,2), Coords(nnj(1),1) - tn(ntj(1),1) -Coords(j,1))

      theta21_2 = atan2(r_i(2)-r_j(2),r_i(1)-r_j(1))&
                  - atan2(Coords(nnj(2),2) - tn(ntj(2),2)-Coords(j,2), Coords(nnj(2),1) - tn(ntj(2),1) -Coords(j,1))

      theta21_3 = atan2(r_i(2)-r_j(2),r_i(1)-r_j(1))&
                  - atan2(Coords(nnj(3),2) - tn(ntj(3),2)-Coords(j,2), Coords(nnj(3),1) - tn(ntj(3),1) -Coords(j,1))


      theta12_1 = atan2(r_i(2)-r_j(2),r_i(1)-r_j(1))&
                  - atan2(Coords(nni(1),2) - tn(nti(1),2)-Coords(i,2), Coords(nni(1),1) - tn(nti(1),1) -Coords(i,1)) + pi

      theta12_2 = atan2(r_i(2)-r_j(2),r_i(1)-r_j(1))&
                  - atan2(Coords(nni(2),2) - tn(nti(2),2)-Coords(i,2), Coords(nni(2),1) - tn(nti(2),1) -Coords(i,1)) + pi

      theta12_3 = atan2(r_i(2)-r_j(2),r_i(1)-r_j(1))&
                  - atan2(Coords(nni(3),2) - tn(nti(3),2)-Coords(i,2), Coords(nni(3),1) - tn(nti(3),1) -Coords(i,1)) + pi

      r=sqrt((r_i(1)-r_j(1))**2 + (r_i(2)-r_j(2))**2)
      V0 = lambda0*exp(-xi0*(r**2))*cos(kappa0*r)
      V3 = lambda3*(r**2.0_dp)*exp(-xi3*(r-y3)**2)
      V6 = lambda6*exp(-xi6*(r-y6)**2)*sin(kappa6*r)

      fTB = V0 + V3*((cos(3.0_dp*theta12_1) + cos(3.0_dp*theta12_2) + cos(3.0_dp*theta12_3))/3.0_dp&
                           + (cos(3.0_dp*theta21_1) + cos(3.0_dp*theta21_2) + cos(3.0_dp*theta21_3))/3.0_dp)&
                      + V6*((cos(6.0_dp*theta12_1) + cos(6.0_dp*theta12_2) + cos(6.0_dp*theta12_3))/3.0_dp&
                           + (cos(6.0_dp*theta21_1) + cos(6.0_dp*theta21_2) + cos(6.0_dp*theta21_3))/3.0_dp)

    else

      r=sqrt((r_i(1)-r_j(1))**2 + (r_i(2)-r_j(2))**2)

      if(r.GT.a0*0.1_dp)then
        fTB = tparallel0*exp(-al0*(r**2))*cos(bet0*r) + tparallel1*(r**2)*exp(-al1*(r-y1)**2)
      else
        fTB=0.0_dp
      endif

    endif

  case default

    write(*,*) 'ERROR: TBFunction must be 1 or 2, got ',TBFunction
    stop 1

  end select

end function fTB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure function fValley(ri,rj,RotMatrix)

  real(dp), intent(in) :: ri(3), rj(3), RotMatrix(2,2)
  real(dp) :: r(2),x,y,phi,rij
  complex(dp) :: fValley

  rij = sqrt((ri(1)-rj(1))**2 + (ri(2)-rj(2))**2)
  if((rij.lt.1.2_dp).and.(rij.gt.0.2_dp))then
      r = matmul(RotMatrix,[ri(1)- rj(1),ri(2) - rj(2)])
      x = r(1)
      y = r(2)

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

  fValley = cmplx(0.0_dp,sign(1.0_dp,cos(3.0_dp*phi)),dp)/sqrt(27.0_dp)
  
     else

  fValley = cmplx(0.0_dp,0.0_dp,dp)

  endif
  
end function fValley

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine PlotBands(Bands,FermiEnergy,MomentaFlattened,numb,numk,Nk,g1,g12)
   
  integer(dp), intent(in)    :: numb,numk,Nk
  integer(dp), intent(in)    :: MomentaFlattened(numk,numk)
  real(dp)   , intent(in)    :: FermiEnergy,g1(2), g12(2), Bands(numb,Nk)
  
  integer(dp) :: i,ivk,ivk1,ivk2,p=0_dp
  real(dp) :: vK1x,vK1y,vK2x,vK2y,vMx,vMy
  real(dp) :: a11,a12,a21,a22,b11,b12,b21,b22,vkx,vky,det
  real(dp) :: aGM, aKG, aMK, aT
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  vK1x=(g1(1)+g12(1))/3.0_dp
  vK1y=(g1(2)+g12(2))/3.0_dp
 
  vK2x=2.0_dp*(g1(1)+g12(1))/3.0_dp
  vK2y=2.0_dp*(g1(2)+g12(2))/3.0_dp
 
  vMx=(g1(1)+g12(1))/2.0_dp
  vMy=(g1(2)+g12(2))/2.0_dp
  
  aKG=sqrt(vK1x**2+vK1y**2)
  aGM=sqrt((vMx-g12(1))**2 + (vMy-g12(2))**2)
  aMK=sqrt((vK1x-vMx)**2 + (vK1y-vMy)**2)

  aT = 3*aKG+aGM+aMK
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  do ivk=0,numk/3-1
    ivk1 = mod(2*numk/3+ivk,numk)
    ivk2 = mod(2*numk/3-2*ivk + p*numk/2,numk)
     write(99,fmt='(F16.8,3X)', advance="no") ivk*aKG/real(numk/3,dp)/aT
     do i=1,numb-1
        write(99,fmt='(F16.8,3X)', advance="no") Bands(i,MomentaFlattened(ivk1+1,ivk2+1))-FermiEnergy
     enddo
     write(99,fmt='(F16.8,3X)') Bands(numb,MomentaFlattened(ivk1+1,ivk2+1))-FermiEnergy
  enddo

  write(99,fmt='(F16.8,3X)', advance="no") aKG/aT
  do i=1,numb-1
     write(99,fmt='(F16.8,3X)', advance="no") Bands(i,MomentaFlattened(1,p*numk/2+1))-FermiEnergy
  enddo
  write(99,fmt='(F16.8,3X)') Bands(numb,MomentaFlattened(1,p*numk/2+1))-FermiEnergy
 
  do ivk=1,numk/3
    ivk1 = mod(ivk,numk)
    ivk2 = mod(numk - 2*ivk + p*numk/2,numk)

    write(99,fmt='(F16.8,3X)', advance="no") aKG/aT + aKG*ivk/real(numk/3,dp)/aT
    do i=1,numb-1
    write(99,fmt='(F16.8,3X)', advance="no") Bands(i,MomentaFlattened(ivk1+1,ivk2+1))-FermiEnergy
    enddo
    write(99,fmt='(F16.8,3X)') Bands(numb,MomentaFlattened(ivk1+1,ivk2+1))-FermiEnergy
 enddo
 
  do ivk=1,numk/6
    ivk1 = mod(numk/3+ivk,numk)
    ivk2 = mod(numk/3+ivk + p*numk/2,numk)

    write(99,fmt='(F16.8,3X)', advance="no") 2*aKG/aT + aMK*ivk/real(numk/6,dp)/aT
    do i=1,numb-1
    write(99,fmt='(F16.8,3X)', advance="no") Bands(i,MomentaFlattened(ivk1+1,ivk2+1))-FermiEnergy
    enddo
    write(99,fmt='(F16.8,3X)') Bands(numb,MomentaFlattened(ivk1+1,ivk2+1))-FermiEnergy
  enddo
 
  do ivk=1,numk/2-1
    ivk1 = mod(numk/2-ivk, numk)
    ivk2 = mod(numk/2+ivk + p*numk/2, numk)

    write(99,fmt='(F16.8,3X)', advance="no") (2*aKG+aMK)/aT + aGM*ivk/real(numk/2,dp)/aT
    do i=1,numb-1
    write(99,fmt='(F16.8,3X)', advance="no") Bands(i,MomentaFlattened(ivk1+1,ivk2+1))-FermiEnergy
    enddo
    write(99,fmt='(F16.8,3X)') Bands(numb,MomentaFlattened(ivk1+1,ivk2+1))-FermiEnergy
  enddo

  do ivk=0,numk/3
    ivk1 = mod(ivk,numk)
    ivk2 = mod(ivk+ p*numk/2,numk)

    write(99,fmt='(F16.8,3X)', advance="no") (2*aKG+aMK+aGM)/aT + ivk*aKG/real(numk/3,dp)/aT
    do i=1,numb-1
    write(99,fmt='(F16.8,3X)', advance="no") Bands(i,MomentaFlattened(ivk1+1,ivk2+1))-FermiEnergy
    enddo
    write(99,fmt='(F16.8,3X)') Bands(numb,MomentaFlattened(ivk1+1,ivk2+1))-FermiEnergy
  enddo
 
  close(99)
  
 end subroutine PlotBands

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dos(Bands,FermiEnergy,MomentaFlattened,numb,numk,Nk,g1,g12)

 integer(dp), intent(in)    :: numb,numk,Nk
 integer(dp), intent(in)    :: MomentaFlattened(numk,numk)
 real(dp)   , intent(in)    :: FermiEnergy,g1(2),g12(2),Bands(numb,Nk)
 
 real(dp) :: OmegaMin,OmegaMax
 integer(dp), parameter :: nspacing=1000  
 integer(dp) :: nomega 
 integer(dp) :: nshift 

 real(dp), allocatable    :: dosCB(:),temp1(:),PivotEnergy(:,:),PivotMomentum(:,:)

 
 integer(dp) :: n,i,ivk,ivk1,ivk2,nbounds(4)
 real(dp) :: Ac,sigma0
 real(dp) :: bandsFull(numb,numk+1,numk+1)
 real(dp) :: vEnergy(numb,numk+1,2),vMomentum(2,numk+1,2)
 
 
  OmegaMin = 0.200_dp
  OmegaMax = -0.200_dp
  
  !  do i=1,Nk
  !  if(maxval(Bands(numb/2-1:numb/2+2,i))-FermiEnergy.GT.OmegaMax)then
  !  OmegaMax=maxval(Bands(numb/2-1:numb/2+2,i))-FermiEnergy
  !  endif
  !  if(minval(Bands(numb/2-1:numb/2+2,i))-FermiEnergy.LT.OmegaMin)then
  !  OmegaMin=minval(Bands(numb/2-1:numb/2+2,i))-FermiEnergy
  !  endif
  !  enddo

  do i=1,Nk
  if(maxval(Bands(:,i))-FermiEnergy.GT.OmegaMax)then
  OmegaMax=maxval(Bands(:,i))-FermiEnergy
  endif
  if(minval(Bands(:,i))-FermiEnergy.LT.OmegaMin)then
  OmegaMin=minval(Bands(:,i))-FermiEnergy
  endif
  enddo
  
  !write(*,*) 'maxmin',OmegaMax,OmegaMin
  
  nomega=int((OmegaMax-OmegaMin)*nspacing)
  nshift=int(-OmegaMin*nspacing)
 
 
 allocate(dosCB(nomega),temp1(nomega))
 
 allocate(pivotEnergy(numb,4),pivotMomentum(2,4))


 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!! Periodic continuation
 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 do ivk1=0,numk-1
 do ivk2=0,numk-1
 bandsFull(:,ivk1+1,ivk2+1)=Bands(:,MomentaFlattened(ivk1+1,ivk2+1))-FermiEnergy
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

  nbounds(1)=1 ! min for i_v
  nbounds(2)=numb-ncb ! max for i_v
  nbounds(3)=numb-ncb+1 ! min for i_c
  nbounds(4)=numb ! max for i_c


  do ivk1=0,numk-1
    do ivk2 = 0,numk
      vEnergy(:,ivk2+1,1)=bandsFull(:,ivk1+1,ivk2+1)
      vEnergy(:,ivk2+1,2)=bandsFull(:,ivk1+2,ivk2+1)
      vMomentum(:,ivk2+1,1)=(ivk1*g1(:)+ivk2*g12(:))/real(numk,dp)
      vMomentum(:,ivk2+1,2)=((ivk1+1)*g1(:)+ivk2*g12(:))/real(numk,dp)
    end do
 
    temp1=0._dp


    !$omp parallel do &
    !$omp private(ivk2,pivotEnergy,pivotMomentum,temp1) &
    !$omp shared(vEnergy,vMomentum,nshift,nomega,nbounds) &
    !$omp reduction(+:dosCB)
    !!$omp reduction(+:dosVB)
    !!$omp reduction(+:dosJoint)  
    do ivk2=0,numk-1

      pivotEnergy=reshape([vEnergy(:,ivk2+1,1),vEnergy(:,ivk2+1,2),vEnergy(:,ivk2+2,1),vEnergy(:,ivk2+2,2)],[numb,4_dp])

      pivotMomentum=reshape([vMomentum(:,ivk2+1,1),vMomentum(:,ivk2+1,2),vMomentum(:,ivk2+2,1),vMomentum(:,ivk2+2,2)],[2_dp,4_dp])

      call dos_S(temp1,pivotMomentum,pivotEnergy,nomega,nspacing,nshift,nbounds)

      !call dosVB_S(ndim,temp2(:),pivotMomentum,pivotEnergy, &
      !nomega,nspacing,nshift,FermiEnergy,nbounds)

      !call dosJoint_S(ndim,temp3(:),pivotMomentum,pivotEnergy, &
      !nomega,nspacing,nshift,nbounds)

      dosCB(:)=dosCB(:)+temp1(:)
      !dosVB(:)=dosVB(:)+temp2(:)
      !dosJoint(:)=dosJoint(:)+temp3(:)


    enddo ! ivk2
    !$omp end parallel do

  enddo ! ivk1

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Ac=4._dp*pi*pi ! normalization: int DOS*sqrt(3)/2*ndim/4= number of bands =nUpper-nLower+1

  ! In units of e^2/hbar=1 and without spin and valley channel
      
  sigma0=1._dp/16._dp

  !do n=nomega,nshift,-1
  
      !write(99,*) -real(n-nshift)/real(nspacing),dosVB(n)/Ac

  !end do
  
  do n=1,nomega
  
    write(99,*) real(n-nshift)/real(nspacing),dosCB(n)/Ac

  end do

  close(99)

 
end subroutine dos

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dos_S(dos,vk,Energies,nomega,nspacing,nshift,nbounds)


 integer(dp), intent(in) :: nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: vk(:,:),Energies(:,:)
 real(dp)   , intent(out) :: dos(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4)
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: add,weight


  dos=0._dp
  do i_c = nbounds(1),nbounds(4)

    energy(:)=Energies(i_c,:)


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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dosVB_S(ndim,dos,vk,Energies,nomega,nspacing,nshift,efermi,nbounds)


 integer(dp), intent(in) :: ndim,nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: efermi
 real(dp)   , intent(in) :: vk(:,:),Energies(:,:)
 real(dp)   , intent(out) :: dos(:)

 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE
 real(dp) :: energy(4)
 real(dp) :: eMax,eMin,dMax,dMin
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: add,weight

 dos=0._dp

 do i_v = nbounds(1),nbounds(2)

   energy(:)=efermi-Energies(i_v,:)


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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dosJoint_S(ndim,dos,vk,Energies,nomega,nspacing,nshift,nbounds)


 integer(dp), intent(in) :: ndim,nomega,nspacing,nshift,nbounds(4)
 real(dp)   , intent(in) :: vk(:,:),Energies(:,:)
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


   energy(:)=Energies(i_c,:)-Energies(i_v,:)

   


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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine n1n2toind(index,n1,n2,nUnitCell_1,nUnitCell_2,numNeighborCells)

    integer(dp), intent(in) :: numNeighborCells
    integer(dp), intent(in) :: n1, n2, nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells)
    integer(dp), intent(inout) :: index

    integer(dp) :: n
    
    index = -1_dp 
    do n=1,numNeighborCells
        if((abs(real(nUnitCell_1(n))-real(n1)).lt.1e-3).and.(abs(real(nUnitCell_2(n))-real(n2)).lt.1e-3))then
            index = n
        endif
    enddo

end subroutine n1n2toind

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dsort_lists(array)

  real(dp) , intent(inout):: array(:,:)
  real(dp), allocatable :: copy(:)

  integer(dp) :: l, ir, i, j
  integer(dp) :: icopy

  l=size(array(:,1))/2+1
  ir=size(array(:,1))

  allocate(copy(size(array(1,:))))

10 continue
  if(l > 1)then
    l=l-1
    copy  =  array(l,:)
  else
    copy  =  array(ir,:)
    array(ir,:)=array(1,:)
    ir=ir-1
    if(ir == 1)then
      array(1,:) = copy
      return
    end if
  end if
  i=l
  j=2*l
20 if(j.le.ir)then
  if(j < ir)then
    if(array(j,1) < array(j+1,1))      j=j+1
  end if
  if(copy(1) < array(j,1))then
    array(i,:) = array(j,:)
    i=j; j=2*j
  else
    j=ir+1
  end if
  goto 20
  end if
  array(i,:) = copy
  goto 10

end subroutine dsort_lists

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dsort_list(array)

  real(dp) , intent(inout):: array(:)
  real(dp) :: copy

  integer(dp) :: l, ir, i, j
  integer(dp) :: icopy

  l=size(array)/2+1
  ir=size(array)


10 continue
  if(l > 1)then
    l=l-1
    copy  =  array(l)
  else
    copy  =  array(ir)
    array(ir)=array(1)
    ir=ir-1
    if(ir == 1)then
      array(1) = copy
      return
    end if
  end if
  i=l
  j=2*l
20 if(j.le.ir)then
  if(j < ir)then
    if(array(j) < array(j+1))      j=j+1
  end if
  if(copy < array(j))then
    array(i) = array(j)
    i=j; j=2*j
  else
    j=ir+1
  end if
  goto 20
  end if
  array(i) = copy
  goto 10

end subroutine dsort_list

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module TightBinding


! Postprocessing compatibility namespace (shares the dev physical settings).
! Compatibility routines grouped by the same responsibilities as FORGE dev.
module ResponseTightBinding
use ResponseSetup
use lapack_routines
implicit none
contains

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


!c!c!cC
!c!c!cC
!c!c!cC


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

end module ResponseTightBinding
