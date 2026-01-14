module TightBinding

use Setup
use lapack_routines

implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine TightBindingHamiltonian(zH,Delta,Coords,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,vk,tn)
 
 integer(dp), intent(in)    :: ndim,numNeighborCells
 integer(dp), intent(in)    :: nUnitCell_1(numNeighborCells),nUnitCell_2(numNeighborCells)
 real(dp)   , intent(in)    :: tn(2,3), vk(2), Coords(ndim,3), Delta
 complex(dp), intent(out) :: zH(ndim,ndim)
 
 integer(dp) :: i, j, n1,n2,icount
 real(dp)    :: rij, Coords_diff
 complex(dp) :: zi_vk_tn,zphase
 
  zH = cmplx(0.0_dp,0.0_dp,dp)

  do i = 1 , ndim
    do j = 1 , i
   
      zphase=exp(cmplx(0.0_dp,-dot_product(vk,Coords(i,1:2)-Coords(j,1:2))*real(fphase,dp),dp))
     
      zH(i,j) = zH(i,j) + ft(Coords(i,:),Coords(j,:))*zphase
      do icount=2,numNeighborCells
        n1=nUnitCell_1(icount)
        n2=nUnitCell_2(icount)
            
        zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
        zH(i,j)=zH(i,j)+ft(Coords(i,:), [Coords(j,1:2)-n1*tn(:,1)-n2*tn(:,2),Coords(j,3)])*exp(-zi_vk_tn)*zphase

                
        zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
        zH(i,j)=zH(i,j)+ft(Coords(i,:), [Coords(j,1:2)+n1*tn(:,1)+n2*tn(:,2),Coords(j,3)])*exp(-zi_vk_tn)*zphase

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
          zV(i,j) = zV(i,j) + ftvalley(ri,rj,transpose(RotMatrix))*zphase
        else if (RotateLayers(nlayer).eq.1)then
          zV(i,j) = zV(i,j) + ftvalley(ri,rj,RotMatrix)*zphase
        endif

        do icount=2,numNeighborCells
          n1=nUnitCell_1(icount)
          n2=nUnitCell_2(icount)
        
          ri =  Coords(i,:)
          rj =  [Coords(j,1:2) - n1*tn(:,1) - n2*tn(:,2), Coords(j,3)]
          
          zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
          
          if (RotateLayers(nlayer).eq.-1)then
            zV(i,j) = zV(i,j) + ftvalley(ri,rj,transpose(RotMatrix))*zphase*exp(-zi_vk_tn)
          else if (RotateLayers(nlayer).eq.1)then
            zV(i,j) = zV(i,j) + ftvalley(ri,rj,RotMatrix)*zphase*exp(-zi_vk_tn)
          endif
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ri =  [Coords(i,1:2) - n1*tn(:,1) - n2*tn(:,2), Coords(i,3)]
          rj =  Coords(j,:) 
          
          zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
          
          if (RotateLayers(nlayer).eq.-1)then
            zV(i,j) = zV(i,j) + ftvalley(ri,rj,transpose(RotMatrix))*zphase*exp(-zi_vk_tn)
          else if (RotateLayers(nlayer).eq.1)then
            zV(i,j) = zV(i,j) + ftvalley(ri,rj,RotMatrix)*zphase*exp(-zi_vk_tn)
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
          zV(i,j) = zV(i,j) + ftvalley(ri,rj,transpose(RotMatrix))*zphase
        else if (RotateLayers(nlayer).eq.1)then
          zV(i,j) = zV(i,j) + ftvalley(ri,rj,RotMatrix)*zphase
        endif

        do icount=2,numNeighborCells
          n1=nUnitCell_1(icount)
          n2=nUnitCell_2(icount)
        
          ri =  Coords(i,:)
          rj =  [Coords(j,1:2) - n1*tn(:,1) - n2*tn(:,2), Coords(j,3)]
          
          zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
          
          if (RotateLayers(nlayer).eq.-1)then
            zV(i,j) = zV(i,j) + ftvalley(ri,rj,transpose(RotMatrix))*zphase*exp(-zi_vk_tn)
          else if (RotateLayers(nlayer).eq.1)then
            zV(i,j) = zV(i,j) + ftvalley(ri,rj,RotMatrix)*zphase*exp(-zi_vk_tn)
          endif         
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ri =  [Coords(i,1:2) - n1*tn(:,1) - n2*tn(:,2), Coords(i,3)]
          rj =  Coords(j,:) 
          
          zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
          
          if (RotateLayers(nlayer).eq.-1)then
            zV(i,j) = zV(i,j) + ftvalley(ri,rj,transpose(RotMatrix))*zphase*exp(-zi_vk_tn)
          else if (RotateLayers(nlayer).eq.1)then
            zV(i,j) = zV(i,j) + ftvalley(ri,rj,RotMatrix)*zphase*exp(-zi_vk_tn)
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

pure function ft(r_i,r_j)

  real(dp), intent(in) :: r_i(3), r_j(3)
  
  real(dp) :: r, cs, ft
  
  r = norm2(r_i-r_j)
  cs = (r_i(3)-r_j(3))/norm2(r_i-r_j)

  ft = 0.0_dp
  if(r.GT.0.001_dp)then
        ft = -2.7_dp*exp((a0-r)/r0)*(1-cs**2) + 0.48_dp*exp((d0-r)/r0)*cs**2
  endif

end function ft

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%

pure function ftvalley(ri,rj,RotMatrix)

  real(dp), intent(in) :: ri(3), rj(3), RotMatrix(2,2)
  real(dp) :: r(2),x,y,phi,rij
  complex(dp) :: ftvalley

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

  ftvalley = cmplx(0.0_dp,sign(1.0_dp,cos(3.0_dp*phi)),dp)/sqrt(27.0_dp)
  
     else

  ftvalley = cmplx(0.0_dp,0.0_dp,dp)

  endif
  
end function ftvalley

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
