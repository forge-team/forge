module HartreeFock

use omp_lib
use Setup
use TightBinding
use lapack_routines

implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine LongRangeInteraction(LongRange,Coords,ndim,t1,t2)

  integer(dp), intent(in) :: ndim
  real(dp), intent(in)    :: t1(2), t2(2),Coords(ndim,3)
  real(dp), intent(inout)    :: LongRange(ndim,ndim)

  integer(dp) :: n,m, n1, n2

  LongRange=0.0_dp

  !$omp parallel do &
  !$omp private(n,m,n1,n2) &
  !$omp shared(Coords,t1,t2) &
  !$omp reduction(+:LongRange)
  do n=1,ndim
    do m=n,ndim

      do n1=-numC,numC
        do n2=-numC,numC

          if(n1+n2.LE.numC.AND.n1+n2.GE.-numC)then
          ! if((norm2(n1*t1+n2*t2)/norm2(t1)).lt.real(numC))then

            LongRange(n,m) = LongRange(n,m) + fv([Coords(n,1:2)-n1*t1-n2*t2,Coords(n,3)], Coords(m,:))

          endif

        enddo
      enddo
    enddo
  enddo
  !$omp end parallel do

  do n=1,ndim-1
    do m=n+1,ndim
      LongRange(m,n)=LongRange(n,m)
    enddo
  enddo

end subroutine LongRangeInteraction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine HamiltonianHartreeFock(zH,Coords,Potential,alpha,Delta,zFock,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,vk,tn)

  integer(dp), intent(in)    :: ndim, numNeighborCells

  integer(dp), intent(in)    :: nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells)
  real(dp)   , intent(in)    :: tn(2,3), vk(2), Coords(ndim,3)
  real(dp)   , intent(in)    :: alpha, Delta, Potential(ndim)
  complex(dp) , intent(in)   :: zFock(ndim,ndim,numNeighborCells)
  complex(dp), intent(out) :: zH(ndim,ndim)

  integer(dp) :: i, j, icount, n1, n2
  real(dp)    :: rij, Coords_diff
  complex(dp) :: zi_vk_tn,zsum,zphase

  zH(:,:) = cmplx(0.0_dp,0.0_dp,dp)

  do i = 1 , ndim
    ! Hartree part
    zH(i,i) = zH(i,i) + Potential(i)
    do j = 1 , i

      ! Kinetic term
      zphase=exp(cmplx(0.0_dp,-dot_product(vk,Coords(i,1:2)-Coords(j,1:2))*real(fphase,dp),dp))

      zH(i,j)=zH(i,j)-alpha*fv(Coords(i,:),Coords(j,:))*conjg(zFock(i,j,1))*zphase
     
      zH(i,j) = zH(i,j)+ft(Coords(i,:), Coords(j,:))*zphase


      do icount=2,numNeighborCells
        n1=nUnitCell_1(icount)
        n2=nUnitCell_2(icount)
              
        zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
        zH(i,j)=zH(i,j)+ft(Coords(i,:), [Coords(j,1:2)-n1*tn(:,1)-n2*tn(:,2),Coords(j,3)])*exp(-zi_vk_tn)*zphase

        
        zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
        zH(i,j)=zH(i,j)+ft(Coords(i,:), [Coords(j,1:2)+n1*tn(:,1)+n2*tn(:,2),Coords(j,3)])*exp(-zi_vk_tn)*zphase

  
      enddo
     
      ! Fock term
      zsum=cmplx(0.0_dp,0.0_dp,dp)
      do icount=2,numNeighborCells
        n1=nUnitCell_1(icount)
        n2=nUnitCell_2(icount)
        
        zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)

        zsum=zsum-alpha*fv(Coords(i,:),[Coords(j,1:2)-n1*tn(:,1)-n2*tn(:,2),Coords(j,3)])*conjg(zFock(i,j,icount))*exp(-zi_vk_tn)
        
        zsum=zsum-alpha*fv(Coords(j,:),[Coords(i,1:2)-n1*tn(:,1)-n2*tn(:,2),Coords(i,3)])*zFock(j,i,icount)*exp( zi_vk_tn)
      
      enddo

      zH(i,j)=zH(i,j)+zsum*zphase
 
    end do
  end do

  ! Layer bias term
  if(nlayers.GT.1)then
  do n1 = 0,nlayers-1
    do i=1,ndim/nlayers
      zH(i+ n1*ndim/nlayers,i+ n1*ndim/nlayers) = zH(i+ n1*ndim/nlayers,i+ n1*ndim/nlayers) + (n1 - real(nlayers-1)/2)*Delta/real(nlayers-1)
    enddo
  enddo
  endif
end subroutine HamiltonianHartreeFock

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Solve_C3(zFockBulk,zFock,Potential,alpha,Bands,zEigenvectors,&
  Coords,MomentaValues,nMomentaComponents,nLower,nUpper,ndim,RotateLayers,numk,Nk,numNeighborCells,nUnitCell_1,nUnitCell_2,tn,g1,g12,nC3pairs)

  integer(dp), intent(in)    :: ndim, numk, Nk, nLower, nUpper,numNeighborCells, nC3pairs(ndim), RotateLayers(nlayers)
  integer(dp), intent(in)    :: nMomentaComponents(Nk,2)
  integer(dp), intent(in)    :: nUnitCell_1(numNeighborCells),nUnitCell_2(numNeighborCells)
  real(dp)   , intent(in)    :: tn(2,3), MomentaValues(numk,numk,2), Coords(ndim,3),Potential(ndim), g1(2), g12(2)
  real(dp)   , intent(in)    :: alpha
  complex(dp), intent(in) :: zFock(ndim,ndim,numNeighborCells)
  real(dp)   , intent(out) :: Bands(nUpper-nLower+1,Nk)
  complex(dp), intent(out) :: zEigenvectors(ndim,nUpper-nLower+1,Nk),zFockBulk(ndim,ndim,numNeighborCells)

  integer(dp) :: ivk1, ivk2, ivk1c3, ivk2c3, icount, icountk, nband, C3_RelatedMomenta(Nk,3)
  integer(dp) :: n1,n2,nlayer,iUnitCell,i, NkMod_C3, NkWith_C3, MomentaMod_C3(2*Nk), MomentaIncl_C3(2*Nk), flg(1)
  real(dp)    :: Energies(ndim), vk(2), vk3(2)
  complex(dp) :: zH_FockBulk1Temp(ndim,ndim)
  complex(dp) :: zP(ndim,1)
  complex(dp) :: zFockBulkTemp(ndim,ndim,numNeighborCells)
  complex(dp) :: zphase(ndim)

  ! Momenta related by C_3
  do icount=1,Nk

    ivk1=nMomentaComponents(icount,1)
    ivk2=nMomentaComponents(icount,2)
    C3_RelatedMomenta(icount,1) = icount    
    
    ivk1c3 = modulo(-ivk1-ivk2,numk)
    ivk2c3 = modulo(ivk1,numk)
    C3_RelatedMomenta(icount,2) = ivk1c3*numk + ivk2c3 +1    
    
    ivk1c3 = modulo(ivk2,numk)
    ivk2c3 = modulo(-ivk1-ivk2,numk)
    C3_RelatedMomenta(icount,3) = ivk1c3*numk + ivk2c3 +1 
        
  enddo

  ! Pick one momentum out of the 3 related by C_3 (or 1 for invariant points)
  MomentaMod_C3(:) = 0_dp
  MomentaIncl_C3(:) = 0_dp
  NkMod_C3 = 0
  NkWith_C3 = 0
  do icount=1,Nk
    flg = findloc(MomentaIncl_C3,icount)
    if (flg(1).eq.0_dp)then
      NkMod_C3 = NkMod_C3+1
      MomentaMod_C3(NkMod_C3) = icount
      
      NkWith_C3 = NkWith_C3+1
      MomentaIncl_C3(NkWith_C3) = icount

      NkWith_C3 = NkWith_C3+1
      MomentaIncl_C3(NkWith_C3) = C3_RelatedMomenta(icount,2)

      NkWith_C3 = NkWith_C3+1
      MomentaIncl_C3(NkWith_C3) = C3_RelatedMomenta(icount,3)
    endif
  enddo

  zFockBulk=cmplx(0.0_dp,0.0_dp,dp)
  !$omp parallel do &
  !$omp private(icount,icountk,vk,vk3,ivk1,ivk2,ivk1c3,ivk2c3,nlayer,nband,zH_FockBulk1Temp,Energies,zP,zphase,iUnitCell,n1,n2,i) &
  !$omp private(zFockBulkTemp) &
  !$omp shared(Coords,Nk,MomentaMod_C3,tn,nLower,nUpper,Bands,zEigenvectors,nUnitCell_1,nUnitCell_2,numNeighborCells,nC3pairs,g1,g12) &
  !$omp shared(nMomentaComponents,MomentaValues,zFock) &
  !$omp reduction(+:zFockBulk)
  do icount=1,NkMod_C3 ! loop over the NkMod_C3 momenta

    icountk = MomentaMod_C3(icount)
    ivk1=nMomentaComponents(icountk,1)
    ivk2=nMomentaComponents(icountk,2)
    vk(:)=MomentaValues(ivk1+1,ivk2+1,:)

    !write(*,*) icount, icountk

    zFockBulkTemp=cmplx(0.0_dp,0.0_dp,dp)

    call HamiltonianHartreeFock(zH_FockBulk1Temp,Coords,Potential,alpha,Delta,zFock,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,vk,tn)

    call diagonalize(zH_FockBulk1Temp,Energies,'V',1_dp,nUpper)

    ! save Bands imposing C_3
    do  nband=nLower,nUpper
      Bands(nband-nLower+1,icountk)=Energies(nband)
      Bands(nband-nLower+1,C3_RelatedMomenta(icountk,2))=Energies(nband)
      Bands(nband-nLower+1,C3_RelatedMomenta(icountk,3))=Energies(nband)
    enddo

    ! Apply C_3 to the central states and save in Eigenvectors
    do nband=nLower,nUpper
      zEigenvectors(:,nband-nLower+1,icountk)=zH_FockBulk1Temp(:,nband)
      
      if ((icountk.ne.1).and.(icountk.ne.(numk/3*(numk+1)+1)).and.(icountk.ne.(2*numk/3*(numk+1)+1)))then
        
        if(fphase.eq.0)then

          ivk1c3 = nMomentaComponents(C3_RelatedMomenta(icountk,2),1)
          ivk2c3 = nMomentaComponents(C3_RelatedMomenta(icountk,2),2)
          zphase(:) = cmplx(1.0_dp,0.0_dp,dp)
          do nlayer=1,nlayers
            if(RotateLayers(nlayer).eq.-1)then
              zphase(ndim/nlayers*(nlayer-1) + ndim/nlayers/2) = exp(cmplx(0.0_dp,-2*pi*real(ivk1c3+ivk2c3,dp)/real(numk,dp),dp))
              zphase(ndim/nlayers*(nlayer-1) + ndim/nlayers) = exp(cmplx(0.0_dp,2*pi*real(ivk1c3+ivk2c3,dp)/real(numk,dp),dp))
            else if (RotateLayers(nlayer).eq.1)then
              zphase(ndim/nlayers*(nlayer-1) + ndim/nlayers/2) = exp(cmplx(0.0_dp,2*pi*real(ivk1c3+ivk2c3,dp)/real(numk,dp),dp))
              zphase(ndim/nlayers*(nlayer-1) + ndim/nlayers) = exp(cmplx(0.0_dp,-2*pi*real(ivk1c3+ivk2c3,dp)/real(numk,dp),dp))         
            endif
          enddo

        elseif(fphase.eq.1)then
          vk3 = real(-ivk1-ivk2 - modulo(-ivk1-ivk2,numk),dp)/real(numk,dp)*g1 +&
            real(ivk1 - modulo(ivk1,numk),dp)/real(numk,dp)*g12
          do i=1,ndim
            zphase(i) = exp(cmplx(0.0_dp,dot_product(vk3,Coords(i,1:2)))) 
          enddo

        endif
        
        zEigenvectors(:,nband-nLower+1,C3_RelatedMomenta(icountk,2)) = zH_FockBulk1Temp(nC3pairs(nC3pairs(:)),nband)*zphase(:)
        
        if(fphase.eq.0)then

          ivk1c3 = nMomentaComponents(C3_RelatedMomenta(icountk,3),1)
          ivk2c3 = nMomentaComponents(C3_RelatedMomenta(icountk,3),2)
          zphase(:) = cmplx(1.0_dp,0.0_dp,dp)
          do nlayer=1,nlayers
            if(RotateLayers(nlayer).eq.-1)then
              zphase(ndim/nlayers*(nlayer-1) + ndim/nlayers/2) = exp(cmplx(0.0_dp,-2*pi*real(ivk2c3,dp)/real(numk,dp),dp))
              zphase(ndim/nlayers*(nlayer-1) + ndim/nlayers) = exp(cmplx(0.0_dp,2*pi*real(ivk2c3,dp)/real(numk,dp),dp))
            else if (RotateLayers(nlayer).eq.1)then
              zphase(ndim/nlayers*(nlayer-1) + ndim/nlayers/2) = exp(cmplx(0.0_dp,2*pi*real(ivk2c3,dp)/real(numk,dp),dp))
              zphase(ndim/nlayers*(nlayer-1) + ndim/nlayers) = exp(cmplx(0.0_dp,-2*pi*real(ivk2c3,dp)/real(numk,dp),dp))         
            endif
          enddo

        elseif(fphase.eq.1)then
          vk3 = real(ivk2 - modulo(ivk2,numk),dp)/real(numk,dp)*g1 +&
            real(-ivk1-ivk2 - modulo(-ivk1-ivk2,numk),dp)/real(numk,dp)*g12
          do i=1,ndim
            zphase(i) = exp(cmplx(0.0_dp,dot_product(vk3,Coords(i,1:2)))) 
          enddo

        endif
        
        zEigenvectors(:,nband-nLower+1,C3_RelatedMomenta(icountk,3)) = zH_FockBulk1Temp(nC3pairs(:),nband)*zphase(:)

      endif
    
    enddo
        
    do i=1,ndim
      zphase(i)=exp(cmplx(0.0_dp,dot_product(vk,Coords(i,1:2))*real(fphase,dp),dp))
    enddo   

    do nband=1,nLower-1
    
      do i=1,ndim
        zP(i,1) = conjg(zH_FockBulk1Temp(i,nband)*zphase(i))
      enddo  
  
      call vector_mul(zP,zFockBulkTemp(:,:,1))
    enddo
      
    do iUnitCell=2,numNeighborCells
      n1=nUnitCell_1(iUnitCell)
      n2=nUnitCell_2(iUnitCell)
      
      zFockBulkTemp(:,:,iUnitCell)=zFockBulkTemp(:,:,1)* &
      exp(cmplx(0.0,-2.0_dp*pi*real(ivk1*n1+ivk2*n1+ivk2*n2,dp)/real(numk,dp),dp))
    enddo
    
    ! Apply C_3 to zFockBulkTemp(:,:,1), obtain the full C_3-tranformed zFockBulkTemp and add the contribution
    ! The points at the corners (ndim/4*m) get a phase related to the unit cell change 
    zH_FockBulk1Temp(:,:) = zFockBulkTemp(nC3pairs(:),nC3pairs(:),1)
    do nlayer=1,nlayers
      if(RotateLayers(nlayer).eq.-1)then

        zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:) =&
              zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:)*exp(-cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 
        zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:) =&
              zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:)*exp(cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp))
        zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2) =&
              zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2)*exp(cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 
        zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers) =&
              zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers)*exp(-cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 

      else if(RotateLayers(nlayer).eq.1)then

        zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:) =&
              zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:)*exp(cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 
        zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:) =&
              zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:)*exp(-cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp))
        zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2) =&
              zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2)*exp(-cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 
        zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers) =&
              zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers)*exp(cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 

      endif
    enddo

    zFockBulkTemp(:,:,1) = zFockBulkTemp(:,:,1) + zH_FockBulk1Temp(:,:)
  
    ivk1c3 = nMomentaComponents(C3_RelatedMomenta(icountk,3),1)
    ivk2c3 = nMomentaComponents(C3_RelatedMomenta(icountk,3),2)
    do iUnitCell=2,numNeighborCells
      n1=nUnitCell_1(iUnitCell)
      n2=nUnitCell_2(iUnitCell)
      
      zFockBulkTemp(:,:,iUnitCell)= zFockBulkTemp(:,:,iUnitCell) + zH_FockBulk1Temp(:,:)* &
      exp(cmplx(0.0,-2.0_dp*pi*real(ivk1c3*n1+ivk2c3*n1+ivk2c3*n2,dp)/real(numk,dp),dp))
      
    enddo
    
    ! Apply C_3^2 to FockBulkTemp(:,:,1) (equivalenlty, C_3 to current FockBulk1Temp(:,:,1)),
    ! obtain the full C_3^2-tranformed FockBulkTemp and add the contribution 
    zH_FockBulk1Temp(:,:) = zH_FockBulk1Temp(nC3pairs(:),nC3pairs(:))
    do nlayer=1,nlayers
      if(RotateLayers(nlayer).eq.-1)then

        zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:) =&
              zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:)*exp(cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 
        zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:) =&
              zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:)*exp(-cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp))
        zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2) =&
              zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2)*exp(-cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 
        zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers) =&
              zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers)*exp(cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 

      else if(RotateLayers(nlayer).eq.1)then

        zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:) =&
              zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:)*exp(-cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 
        zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:) =&
              zH_FockBulk1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:)*exp(cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp))
        zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2) =&
              zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2)*exp(cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 
        zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers) =&
              zH_FockBulk1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers)*exp(-cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 

      endif
    enddo
    
    zFockBulkTemp(:,:,1) = zFockBulkTemp(:,:,1) + zH_FockBulk1Temp(:,:)
    
    ivk1c3 = nMomentaComponents(C3_RelatedMomenta(icountk,2),1)
    ivk2c3 = nMomentaComponents(C3_RelatedMomenta(icountk,2),2)
    do iUnitCell=2,numNeighborCells
      n1=nUnitCell_1(iUnitCell)
      n2=nUnitCell_2(iUnitCell)
      
      zFockBulkTemp(:,:,iUnitCell)= zFockBulkTemp(:,:,iUnitCell) + zH_FockBulk1Temp(:,:)* &
      exp(cmplx(0.0,-2.0_dp*pi*real(ivk1c3*n1+ivk2c3*n1+ivk2c3*n2,dp)/real(numk,dp),dp))
      
    enddo

    ! Add to zFockBulk. If the momentum is C_3-invariant we divide by 3 to conserve the number of states
    if ((icountk.ne.1).and.(icountk.ne.(numk/3*(numk+1)+1)).and.(icountk.ne.(2*numk/3*(numk+1)+1)))then
      zFockBulk=zFockBulk + zFockBulkTemp
    else
      zFockBulk=zFockBulk + zFockBulkTemp/3.0_dp
    endif

  enddo
  !$omp end parallel do

  !zFockBulk=zFockBulk/cmplx(real(numk*numk,dp),0.0_dp,dp)

end subroutine Solve_C3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Solve(zFockBulk,zFock,Potential,alpha,Bands,zEigenvectors,&
  Coords,MomentaValues,nMomentaComponents,nLower,nUpper,ndim,numk,Nk,numNeighborCells,nUnitCell_1,nUnitCell_2,tn)

  integer(dp), intent(in)    :: ndim, numk, Nk, nLower, nUpper,numNeighborCells
  integer(dp), intent(in)    :: nMomentaComponents(Nk,2)
  integer(dp), intent(in)    :: nUnitCell_1(numNeighborCells),nUnitCell_2(numNeighborCells)
  real(dp)   , intent(in)    :: tn(2,3), MomentaValues(numk,numk,2), Coords(ndim,3),Potential(ndim),alpha
  complex(dp), intent(in) :: zFock(ndim,ndim,numNeighborCells)
  real(dp)   , intent(out) :: Bands(nUpper-nLower+1,Nk)
  complex(dp), intent(out) :: zEigenvectors(ndim,nUpper-nLower+1,Nk),zFockBulk(ndim,ndim,numNeighborCells)

  integer(dp) :: ivk1, ivk2, icount, nband
  integer(dp) :: n1,n2,iUnitCell,i
  real(dp)    :: Energies(ndim), vk(2)
  complex(dp) :: H(ndim,ndim)
  complex(dp) :: zP(ndim,1)
  complex(dp) :: zFockBulkTemp(ndim,ndim,numNeighborCells)
  complex(dp) :: zphase(ndim)


  zFockBulk=cmplx(0.0_dp,0.0_dp,dp)

  !$omp parallel do &
  !$omp private(icount,vk,ivk1,ivk2,nband,H,Energies,zP,zphase,iUnitCell,n1,n2,i) &
  !$omp private(zFockBulkTemp) &
  !$omp shared(Coords,Nk,tn,nLower,nUpper,Bands,zEigenvectors,nUnitCell_1,nUnitCell_2,numNeighborCells) &
  !$omp shared(nMomentaComponents,MomentaValues,zFock) &
  !$omp reduction(+:zFockBulk)
  do icount=1,Nk

      !write(*,*) icount

      zFockBulkTemp=cmplx(0.0_dp,0.0_dp,dp)

      ivk1=nMomentaComponents(icount,1)
      ivk2=nMomentaComponents(icount,2)
      vk(:)=MomentaValues(ivk1+1,ivk2+1,:)

      call HamiltonianHartreeFock(H,Coords,Potential,alpha,Delta,zFock,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,vk,tn)

      call diagonalize(H,Energies,'V',1_dp,nUpper)

      do  nband=nLower,nUpper
        Bands(nband-nLower+1,icount)=Energies(nband)
      enddo

      do nband=nLower,nUpper
        zEigenvectors(:,nband-nLower+1,icount)=H(:,nband)
      enddo
          
      do i=1,ndim
        zphase(i)=exp(cmplx(0.0_dp,dot_product(vk,Coords(i,1:2))*real(fphase,dp),dp))
      enddo   
    
      do nband=1,nLower-1
              
        do i=1,ndim
          zP(i,1) = conjg(H(i,nband)*zphase(i))
        enddo  
        
        call vector_mul(zP,zFockBulkTemp(:,:,1))
      
      enddo
      
      do iUnitCell=2,numNeighborCells
        n1=nUnitCell_1(iUnitCell)
        n2=nUnitCell_2(iUnitCell)
          
        zFockBulkTemp(:,:,iUnitCell)=zFockBulkTemp(:,:,1)* &
        exp(cmplx(0.0,-2.0_dp*pi*real(ivk1*n1+ivk2*n1+ivk2*n2,dp)/real(numk,dp),dp))
      
      enddo
      
      zFockBulk=zFockBulk+zFockBulkTemp

  enddo
  !$omp end parallel do

  !zFockBulk=zFockBulk/cmplx(real(numk*numk,dp),0.0_dp,dp)


end subroutine Solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GetFock(zFock,zSortedEigenvectors,ndim,numb,numk,DegFactor,OccStates,PartOccStates,Coords,MomentaValues,numNeighborCells,nUnitCell_1,nUnitCell_2,nSortedMomenta)

  integer(dp), intent(in)    :: ndim,numb,numk,OccStates,PartOccStates,numNeighborCells
  integer(dp), intent(in)    :: nUnitCell_1(numNeighborCells),nUnitCell_2(numNeighborCells),nSortedMomenta(numb*numk*numk,2)
  real(dp), intent(in)        :: Coords(ndim,3),MomentaValues(numk,numk,2),DegFactor
  complex(dp), intent(in) :: zSortedEigenvectors(ndim,numb*numk*numk)
  complex(dp) , intent(inout)   :: zFock(ndim,ndim,numNeighborCells)

  integer(dp) :: n1,n2,ivk1,ivk2, i,nband,iUnitCell
  real(dp) :: vk(2)
  complex(dp) :: zP(ndim,1),zdeg

  !zFock=cmplx(0.0_dp,0.0_dp,dp)
  
  !$omp parallel do & 
  !$omp private(nband, ivk1,ivk2,vk,zP) &
  !$omp private(i,iUnitCell,n1,n2) &
  !$omp shared(OccStates,PartOccStates,nSortedMomenta,MomentaValues,Coords,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,numk,zSortedEigenvectors) &
  !$omp reduction(+:zFock)
  do nband=1,OccStates-PartOccStates
        
    ivk1=nSortedMomenta(nband,1)
    ivk2=nSortedMomenta(nband,2)
    vk(:)=MomentaValues(ivk1+1,ivk2+1,:)
        
   
    do i=1,ndim
      zP(i,1) = conjg(zSortedEigenvectors(i,nband)*exp(cmplx(0.0_dp,dot_product(vk,Coords(i,1:2))*real(fphase,dp),dp)))
    enddo

    do iUnitCell=1,numNeighborCells
      n1=nUnitCell_1(iUnitCell)
      n2=nUnitCell_2(iUnitCell)

      call vector_mul(zP,zFock(:,:,iUnitCell),exp(cmplx(0.0_dp,-2.0_dp*pi*real(ivk1*n1+ivk2*n1+ivk2*n2,dp)/real(numk,dp),dp)))
      
    enddo
 
  enddo
  !$omp end parallel do

  ! Select all degenerate states and choose the average to comply with the correct symmetry        
  zdeg=cmplx(DegFactor,dp)

  !$omp parallel do & 
  !$omp private(nband, ivk1,ivk2,vk,zP) &
  !$omp private(i,iUnitCell,n1,n2) &
  !$omp shared(OccStates,PartOccStates,nSortedMomenta,MomentaValues,Coords,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,numk,zSortedEigenvectors,zdeg) &
  !$omp reduction(+:zFock)
  do nband=OccStates-PartOccStates+1,OccStates

    ivk1=nSortedMomenta(nband,1)
    ivk2=nSortedMomenta(nband,2)
    vk(:)=MomentaValues(ivk1+1,ivk2+1,:)

    do i=1,ndim
      zP(i,1) = conjg(zSortedEigenvectors(i,nband)*exp(cmplx(0.0_dp,dot_product(vk,Coords(i,1:2))*real(fphase,dp),dp)))
    enddo

    do iUnitCell=1,numNeighborCells
      n1=nUnitCell_1(iUnitCell)
      n2=nUnitCell_2(iUnitCell)

      call vector_mul(zP,zFock(:,:,iUnitCell),&
      exp(cmplx(0.0_dp,-2.0_dp*pi*real(ivk1*n1+ivk2*n1+ivk2*n2,dp)/real(numk,dp),dp))*zdeg)

    enddo
  
  enddo
  !$omp end parallel do

  !zFock=zFock/cmplx(real(numk*numk,dp),0.0_dp,dp)

end subroutine GetFock

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GetFock_C3(zFock,zSortedEigenvectors,ndim,RotateLayers,numb,numk,DegFactor,OccStates,PartOccStates,Coords,MomentaValues,numNeighborCells,nUnitCell_1,nUnitCell_2,nSortedMomenta,nC3pairs)

  integer(dp), intent(in)    :: ndim,RotateLayers(nlayers),numb,numk,OccStates,PartOccStates,numNeighborCells,nC3pairs(ndim)
  integer(dp), intent(in)    :: nUnitCell_1(numNeighborCells),nUnitCell_2(numNeighborCells),nSortedMomenta(numb*numk*numk,2)
  real(dp), intent(in)        :: Coords(ndim,3),MomentaValues(numk,numk,2),DegFactor
  complex(dp), intent(in) :: zSortedEigenvectors(ndim,numb*numk*numk)
  complex(dp) , intent(inout)   :: zFock(ndim,ndim,numNeighborCells)

  integer(dp) :: n1,n2,nlayer,ivk1,ivk2, i,nband,iUnitCell,icount
  real(dp) :: vk(2)
  complex(dp) :: zP(ndim,1),zdeg
  complex(dp) :: zFock1Temp(ndim,ndim), zFockTemp(ndim,ndim,numNeighborCells)

  !zFock=cmplx(0.0_dp,0.0_dp,dp)

  !$omp parallel do & 
  !$omp private(nband, ivk1,ivk2,icount,vk,zP) &
  !$omp private(i,iUnitCell,n1,n2) &
  !$omp shared(OccStates,PartOccStates,nSortedMomenta,MomentaValues,Coords,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,numk,zSortedEigenvectors) &
  !$omp reduction(+:zFock)
  do nband=1,OccStates-PartOccStates

    ivk1=nSortedMomenta(nband,1)
    ivk2=nSortedMomenta(nband,2)
    vk(:)=MomentaValues(ivk1+1,ivk2+1,:)

    icount= ivk1*numk + ivk2 + 1

    ! C_3 symmetry is respected on the contributions from C_3-non invariant points, after using Solve_C3
    if((icount.ne.1).and.(icount.ne.(numk/3*(numk+1)+1)).and.(icount.ne.(2*numk/3*(numk+1)+1)))then

      do i=1,ndim
        zP(i,1) = conjg(zSortedEigenvectors(i,nband)*exp(cmplx(0.0_dp,dot_product(vk,Coords(i,1:2))*real(fphase,dp),dp)))
      enddo

      do iUnitCell=1,numNeighborCells
        n1=nUnitCell_1(iUnitCell)
        n2=nUnitCell_2(iUnitCell)

        call vector_mul(zP,zFock(:,:,iUnitCell),exp(cmplx(0.0_dp,-2.0_dp*pi*real(ivk1*n1+ivk2*n1+ivk2*n2,dp)/real(numk,dp),dp)))
        
      enddo
    
    endif

    
  enddo
  !$omp end parallel do

  ! Select all degenerate states and choose the average to comply with the correct symmetry        
  zdeg=cmplx(DegFactor,dp)
  
  !$omp parallel do & 
  !$omp private(nband, ivk1,ivk2,icount,vk,zP) &
  !$omp private(i,iUnitCell,n1,n2) &
  !$omp shared(OccStates,PartOccStates,nSortedMomenta,MomentaValues,Coords,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,numk,zSortedEigenvectors,zdeg) &
  !$omp reduction(+:zFock)
  do nband=OccStates-PartOccStates+1,OccStates

    ivk1=nSortedMomenta(nband,1)
    ivk2=nSortedMomenta(nband,2)
    vk(:)=MomentaValues(ivk1+1,ivk2+1,:)
    
    icount= ivk1*numk + ivk2 + 1

    if((icount.ne.1).and.(icount.ne.(numk/3*(numk+1)+1)).and.(icount.ne.(2*numk/3*(numk+1)+1)))then
      
      do i=1,ndim
        zP(i,1) = conjg(zSortedEigenvectors(i,nband)*exp(cmplx(0.0_dp,dot_product(vk,Coords(i,1:2))*real(fphase,dp),dp)))
      enddo

      do iUnitCell=1,numNeighborCells
        n1=nUnitCell_1(iUnitCell)
        n2=nUnitCell_2(iUnitCell)

        call vector_mul(zP,zFock(:,:,iUnitCell),&
        exp(cmplx(0.0_dp,-2.0_dp*pi*real(ivk1*n1+ivk2*n1+ivk2*n2,dp)/real(numk,dp),dp))*zdeg)

      enddo

    endif
  
  enddo
  !$omp end parallel do

  !$omp parallel do & 
  !$omp private(nband,nlayer, ivk1,ivk2,icount,vk,zP,zFockTemp,zFock1Temp) &
  !$omp private(i,iUnitCell,n1,n2) &
  !$omp shared(OccStates,PartOccStates,nSortedMomenta,MomentaValues,Coords,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,numk,zSortedEigenvectors,nC3pairs) &
  !$omp reduction(+:zFock)
  do nband=1,OccStates-PartOccStates

    ivk1=nSortedMomenta(nband,1)
    ivk2=nSortedMomenta(nband,2)
    vk(:)=MomentaValues(ivk1+1,ivk2+1,:)

    icount= ivk1*numk + ivk2 + 1

    ! We enforce C_3 
    if((icount.eq.1).or.(icount.eq.(numk/3*(numk+1)+1)).or.(icount.eq.(2*numk/3*(numk+1)+1)))then

      zFockTemp(:,:,:) = cmplx(0.0_dp,0.0_dp,dp)
      zFock1Temp(:,:) = cmplx(0.0_dp,0.0_dp,dp)

      do i=1,ndim
        zP(i,1) = conjg(zSortedEigenvectors(i,nband)*exp(cmplx(0.0_dp,dot_product(vk,Coords(i,1:2))*real(fphase,dp),dp)))
      enddo 
            
      call vector_mul(zP,zFockTemp(:,:,1))
            
      do iUnitCell=2,numNeighborCells
        n1=nUnitCell_1(iUnitCell)
        n2=nUnitCell_2(iUnitCell)
        
        zFockTemp(:,:,iUnitCell)=zFockTemp(:,:,1)* &
        exp(cmplx(0.0,-2.0_dp*pi*real(ivk1*n1+ivk2*n1+ivk2*n2,dp)/real(numk,dp),dp))
      enddo
          
      zFock1Temp(:,:) = zFockTemp(nC3pairs(:),nC3pairs(:),1)
      do nlayer=1,nlayers
        if(RotateLayers(nlayer).eq.-1)then

          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:)*exp(-cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 
          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:)*exp(cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp))
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2)*exp(cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers)*exp(-cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 

        else if(RotateLayers(nlayer).eq.1)then

          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:)*exp(cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 
          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:)*exp(-cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp))
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2)*exp(-cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers)*exp(cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 

        endif
      enddo
      
      zFockTemp(:,:,1) = zFockTemp(:,:,1) + zFock1Temp(:,:)
      do iUnitCell=2,numNeighborCells
        n1=nUnitCell_1(iUnitCell)
        n2=nUnitCell_2(iUnitCell)
        
        zFockTemp(:,:,iUnitCell)= zFockTemp(:,:,iUnitCell) + zFock1Temp(:,:)* &
        exp(cmplx(0.0,-2.0_dp*pi*real(ivk1*n1+ivk2*n1+ivk2*n2,dp)/real(numk,dp),dp))
        
      enddo
   
      zFock1Temp(:,:) = zFock1Temp(nC3pairs(:),nC3pairs(:))
      do nlayer=1,nlayers
        if(RotateLayers(nlayer).eq.-1)then

          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:)*exp(cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 
          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:)*exp(-cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp))
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2)*exp(-cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers)*exp(cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 

        else if(RotateLayers(nlayer).eq.1)then

          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:)*exp(-cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 
          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:)*exp(cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp))
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2)*exp(cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers)*exp(-cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 

        endif
      enddo
      
      zFockTemp(:,:,1) = zFockTemp(:,:,1) + zFock1Temp(:,:)
      do iUnitCell=2,numNeighborCells
        n1=nUnitCell_1(iUnitCell)
        n2=nUnitCell_2(iUnitCell)
        
        zFockTemp(:,:,iUnitCell)= zFockTemp(:,:,iUnitCell) + zFock1Temp(:,:)* &
        exp(cmplx(0.0,-2.0_dp*pi*real(ivk1*n1+ivk2*n1+ivk2*n2,dp)/real(numk,dp),dp))
    
      enddo

      zFock(:,:,:) = zFock(:,:,:) + zFockTemp(:,:,:)/cmplx(3.0_dp,0.0_dp,dp)

    endif

  enddo
  !$omp end parallel do

    
  !$omp parallel do & 
  !$omp private(nband,nlayer,ivk1,ivk2,icount,vk,zP,zFockTemp,zFock1Temp) &
  !$omp private(i,iUnitCell,n1,n2) &
  !$omp shared(OccStates,PartOccStates,nSortedMomenta,MomentaValues,Coords,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,numk,zSortedEigenvectors,zdeg,nC3pairs) &
  !$omp reduction(+:zFock)
  do nband=OccStates-PartOccStates+1,OccStates

    ivk1=nSortedMomenta(nband,1)
    ivk2=nSortedMomenta(nband,2)
    vk(:)=MomentaValues(ivk1+1,ivk2+1,:)
    
    icount= ivk1*numk + ivk2 + 1

    ! We enforce C_3 symmetry on the C_3-invariant points, similarly done in Solve_C3
    if((icount.eq.1).or.(icount.eq.(numk/3*(numk+1)+1)).or.(icount.eq.(2*numk/3*(numk+1)+1)))then

      zFockTemp(:,:,:) = cmplx(0.0_dp,0.0_dp,dp)
      zFock1Temp(:,:) = cmplx(0.0_dp,0.0_dp,dp)

      do i=1,ndim
        zP(i,1) = conjg(zSortedEigenvectors(i,nband)*exp(cmplx(0.0_dp,dot_product(vk,Coords(i,1:2))*real(fphase,dp),dp)))
      enddo 
            
      call vector_mul(zP,zFockTemp(:,:,1),zdeg)
            
      do iUnitCell=2,numNeighborCells
        n1=nUnitCell_1(iUnitCell)
        n2=nUnitCell_2(iUnitCell)
        
        zFockTemp(:,:,iUnitCell)=zFockTemp(:,:,1)* &
        exp(cmplx(0.0,-2.0_dp*pi*real(ivk1*n1+ivk2*n1+ivk2*n2,dp)/real(numk,dp),dp))
      enddo
          
      zFock1Temp(:,:) = zFockTemp(nC3pairs(:),nC3pairs(:),1)
      do nlayer=1,nlayers
        if(RotateLayers(nlayer).eq.-1)then

          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:)*exp(-cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 
          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:)*exp(cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp))
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2)*exp(cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers)*exp(-cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 

        else if(RotateLayers(nlayer).eq.1)then

          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:)*exp(cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 
          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:)*exp(-cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp))
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2)*exp(-cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers)*exp(cmplx(0.0_dp,2*pi*real(ivk1+ivk2)/real(numk),dp)) 

        endif
      enddo
      
      zFockTemp(:,:,1) = zFockTemp(:,:,1) + zFock1Temp(:,:)
      do iUnitCell=2,numNeighborCells
        n1=nUnitCell_1(iUnitCell)
        n2=nUnitCell_2(iUnitCell)
        
        zFockTemp(:,:,iUnitCell)= zFockTemp(:,:,iUnitCell) + zFock1Temp(:,:)* &
        exp(cmplx(0.0,-2.0_dp*pi*real(ivk1*n1+ivk2*n1+ivk2*n2,dp)/real(numk,dp),dp))
        
      enddo
      
      zFock1Temp(:,:) = zFock1Temp(nC3pairs(:),nC3pairs(:))
      do nlayer=1,nlayers
        if(RotateLayers(nlayer).eq.-1)then

          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:)*exp(cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 
          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:)*exp(-cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp))
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2)*exp(-cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers)*exp(cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 

        else if(RotateLayers(nlayer).eq.1)then

          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers/2,:)*exp(-cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 
          zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:) =&
                zFock1Temp(ndim/nlayers*(nlayer-1) + ndim/nlayers,:)*exp(cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp))
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers/2)*exp(cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 
          zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers) =&
                zFock1Temp(:,ndim/nlayers*(nlayer-1) + ndim/nlayers)*exp(-cmplx(0.0_dp,2*pi*real(ivk1)/real(numk),dp)) 

        endif
      enddo
      
      zFockTemp(:,:,1) = zFockTemp(:,:,1) + zFock1Temp(:,:)
      do iUnitCell=2,numNeighborCells
        n1=nUnitCell_1(iUnitCell)
        n2=nUnitCell_2(iUnitCell)
        
        zFockTemp(:,:,iUnitCell)= zFockTemp(:,:,iUnitCell) + zFock1Temp(:,:)* &
        exp(cmplx(0.0,-2.0_dp*pi*real(ivk1*n1+ivk2*n1+ivk2*n2,dp)/real(numk,dp),dp))
        
      enddo
      
      ! Divide by 3 to maintain number of states
      zFock(:,:,:) = zFock(:,:,:) + zFockTemp(:,:,:)/cmplx(3.0_dp,0.0_dp,dp)

    endif
  
  enddo
  !$omp end parallel do

  !zFock=zFock/cmplx(real(numk*numk,dp),0.0_dp,dp)

end subroutine GetFock_C3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SortEnergies_1Spin(ndim,numk,numb,Nk,nMomentaComponents,zSortedEigenvectors,SortedEnergies_1Spin,nSortedMomenta,zEigenvectors,Bands)

    integer(dp), intent(in)    :: numb, numk, Nk, ndim
    integer(dp), intent(in)    :: nMomentaComponents(Nk,2)
    
    real(dp)   , intent(in)    :: Bands(numb,Nk)
    complex(dp), intent(in)    :: zEigenvectors(ndim,numb,Nk)
    
    integer(dp), intent(inout) :: nSortedMomenta(numb*numk*numk,2)
    real(dp) , intent(inout)   :: SortedEnergies_1Spin(numb*numk*numk)
    complex(dp), intent(inout) :: zSortedEigenvectors(ndim,numb*numk*numk)

    
    integer(dp) :: i, indev, numNeighborCells, icount, nbands
    real(dp) :: sum, energies(numk*numk*numb,2)
   
    
    do icount=1,numk*numk*numb
        indev = (icount-1_dp)/(numk*numk) + 1_dp
        energies(icount,1) = Bands(indev,icount - (indev-1_dp)*numk*numk)
        energies(icount,2) = real(icount,dp)
    enddo
    
    SortedEnergies_1Spin(:) = 0.0_dp
    zSortedEigenvectors(:,:) = cmplx(0.0_dp,0.0_dp,dp)
    nSortedMomenta(:,:) = 0_dp
    call dsort_lists(energies)
    SortedEnergies_1Spin(:) =  energies(:,1)
    !!$omp parallel do &
    !!$omp private(icount,nbands,nsort) &
    !!$omp shared(ncount1,WindowSort_2Spins,SortedEnergies_1Spin,zSortedEigenvectors,Bands,zEigenvectors,nSortedMomenta,nMomentaComponents)
    do icount=1,numk*numk*numb
        numNeighborCells = nint(energies(icount,2),dp)
        indev = (numNeighborCells-1_dp)/(numk*numk) + 1_dp
        zSortedEigenvectors(:,icount)=zEigenvectors(:,indev,numNeighborCells - (indev-1)*numk*numk)
        nSortedMomenta(icount,1) = nMomentaComponents(numNeighborCells - (indev-1_dp)*numk*numk,1)
        nSortedMomenta(icount,2) = nMomentaComponents(numNeighborCells - (indev-1_dp)*numk*numk,2)
    enddo
    !!$omp end parallel do
    
end subroutine SortEnergies_1Spin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SortEnergies_2Spins(numk,numb,WindowSort_2Spins,SortedEnergies_1Spin,SortedEnergies_2Spins)

  integer(dp), intent(in)    :: numb, numk
  integer(dp), intent(in)    :: WindowSort_2Spins(2,2)

  real(dp) , intent(in)   :: SortedEnergies_1Spin(numb*numk*numk,2)
  real(dp), intent(inout) :: SortedEnergies_2Spins(2*numb*numk*numk,2)

  integer(dp) :: ncountBoth,ncount1,ncount2, nbands

  ncountBoth=0
  ncount1=1
  ncount2=1
  do nbands=1,(WindowSort_2Spins(1,1)-1)*numk*numk
    ncountBoth=ncountBoth+1
    SortedEnergies_2Spins(ncountBoth,1)=SortedEnergies_1Spin(ncount1,1)
    SortedEnergies_2Spins(ncountBoth,2)=1
    ncount1=ncount1+1
  enddo
  do nbands=1,(WindowSort_2Spins(1,2)-1)*numk*numk
    ncountBoth=ncountBoth+1
    SortedEnergies_2Spins(ncountBoth,1)=SortedEnergies_1Spin(ncount1,2)
    SortedEnergies_2Spins(ncountBoth,2)=2
    ncount2=ncount2+1
  enddo

  do while(ncount1.LE.numk*numk*(WindowSort_2Spins(2,1)-1).AND.(ncount2.LE.numk*numk*(WindowSort_2Spins(2,2)-1)))
    if(SortedEnergies_1Spin(ncount1,1).LE.SortedEnergies_1Spin(ncount2,2))then
      ncountBoth=ncountBoth+1
      SortedEnergies_2Spins(ncountBoth,1)=SortedEnergies_1Spin(ncount1,1)
      SortedEnergies_2Spins(ncountBoth,2)=1
      ncount1=ncount1+1
      do while(SortedEnergies_1Spin(ncount1,1).LE.SortedEnergies_1Spin(ncount2,2).AND.(ncount1.LE.numk*numk*numb))
        ncountBoth=ncountBoth+1
        SortedEnergies_2Spins(ncountBoth,1)=SortedEnergies_1Spin(ncount1,1)
        SortedEnergies_2Spins(ncountBoth,2)=1
        ncount1=ncount1+1
      enddo
      ncountBoth=ncountBoth+1
      SortedEnergies_2Spins(ncountBoth,1)=SortedEnergies_1Spin(ncount2,2)
      SortedEnergies_2Spins(ncountBoth,2)=2
      ncount2=ncount2+1
    else
      ncountBoth=ncountBoth+1
      SortedEnergies_2Spins(ncountBoth,1)=SortedEnergies_1Spin(ncount2,2)
      SortedEnergies_2Spins(ncountBoth,2)=2
      ncount2=ncount2+1
      do while(SortedEnergies_1Spin(ncount2,2).LE.SortedEnergies_1Spin(ncount1,1).AND.(ncount2.LE.numk*numk*numb))
        ncountBoth=ncountBoth+1
        SortedEnergies_2Spins(ncountBoth,1)=SortedEnergies_1Spin(ncount2,2)
        SortedEnergies_2Spins(ncountBoth,2)=2
        ncount2=ncount2+1
      enddo
      ncountBoth=ncountBoth+1
      SortedEnergies_2Spins(ncountBoth,1)=SortedEnergies_1Spin(ncount1,1)
      SortedEnergies_2Spins(ncountBoth,2)=1
      ncount1=ncount1+1
    endif
  enddo

  do while(ncount1.LE.numk*numk*numb)
    ncountBoth=ncountBoth+1
    SortedEnergies_2Spins(ncountBoth,1)=SortedEnergies_1Spin(ncount1,1)
    SortedEnergies_2Spins(ncountBoth,2)=1
    ncount1=ncount1+1
  enddo

  do while(ncount2.LE.numk*numk*numb)
    ncountBoth=ncountBoth+1
    SortedEnergies_2Spins(ncountBoth,1)=SortedEnergies_1Spin(ncount2,2)
    SortedEnergies_2Spins(ncountBoth,2)=2
    ncount2=ncount2+1
  enddo

  return
end subroutine SortEnergies_2Spins

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c!c!cC
!c!c!cC
!c!c!cC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GetKineticEnergy(fsum,Coords,zFock,Delta,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,tn)

  integer(dp), intent(in)    :: ndim, numNeighborCells

  integer(dp), intent(in)    :: nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells)
  real(dp)   , intent(in)    :: tn(2,3), Coords(ndim,3), Delta
  complex(dp) , intent(in)   :: zFock(ndim,ndim,numNeighborCells)
  real(dp), intent(out) :: fsum

  integer(dp) :: i, j, icount, n1, n2
  complex(dp) :: zsum

  zsum=cmplx(0.0_dp,0.0_dp,dp)

  !$omp parallel do &
  !$omp private(i,j,icount,n1,n2) &
  !$omp shared(ndim,zFock,numNeighborCells,Coords,tn,nUnitCell_1,nUnitCell_2) &
  !$omp reduction(+:zsum) 
  do i = 1,ndim
   do j = 1,ndim
   
     icount=1
      
     zsum=zsum+ft(Coords(i,:),Coords(j,:))*zFock(i,j,1)

      do icount=2,numNeighborCells
        n1=nUnitCell_1(icount)
        n2=nUnitCell_2(icount)
        
        zsum=zsum+ft(Coords(i,:), [Coords(j,1:2)-n1*tn(:,1)-n2*tn(:,2),Coords(j,3)])*zFock(i,j,icount)

        zsum=zsum+ft(Coords(i,:), [Coords(j,1:2)+n1*tn(:,1)+n2*tn(:,2),Coords(j,3)])*conjg(zFock(j,i,icount))

      enddo

   end do
  end do
  !$omp end parallel do

  do n1 = 0,nlayers-1
    do i=1,ndim/nlayers
       zsum = zsum + zFock(i + n1*ndim/nlayers,i + n1*ndim/nlayers,1)*(n1 - real(nlayers-1)/2)*Delta/real(nlayers-1)
    enddo
  enddo

  fsum=real(zsum,dp)

end subroutine GetKineticEnergy


!c!c!cC
!c!c!cC
!c!c!cC

subroutine GetHartreeHubbardEnergy(fsumH,fsumU,LongRange,dens1,dens2,ndim)

  integer(dp), intent(in)    :: ndim
  real(dp), intent(in)    :: LongRange(ndim,ndim)
  real(dp), intent(in)  :: dens1(ndim),dens2(ndim)
  real(dp), intent(out) :: fsumH,fsumU

  integer(dp) :: n,m

  fsumH=0.0_dp
  !$omp parallel do &
  !$omp private(n,m) &
  !$omp shared(dens1,dens2,LongRange)&
  !$omp reduction(+:fsumH)
  do n=1,ndim
    do m=1,ndim
      fsumH = fsumH + 0.5_dp*(dens1(m)+dens2(m))*(dens1(n)+dens2(n))*LongRange(n,m)
    enddo
  enddo
  !$omp end parallel do

  fsumU=0.0_dp
  do n=1,ndim
    fsumU = fsumU + dens1(n)*dens2(n)
  enddo

end subroutine GetHartreeHubbardEnergy

!c!c!cC
!c!c!cC
!c!c!cC

subroutine GetFockEnergy(fsum,Coords,zFock,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,tn)

    integer(dp), intent(in)    :: ndim, numNeighborCells

    integer(dp), intent(in)    :: nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells)
    real(dp)   , intent(in)    :: tn(2,3), Coords(ndim,3)
    complex(dp), intent(in)  :: zFock(ndim,ndim,numNeighborCells)
    real(dp), intent(out) :: fsum

    integer(dp) :: i, j, nt, icount, n1, n2
    real(dp)    :: rij
    complex(dp) :: zsum


   zsum=cmplx(0.0_dp,0.0_dp,dp)

   !$omp parallel do &
   !$omp private(i,j,icount,rij,n1,n2) &
   !$omp shared(ndim,zFock,numNeighborCells,Coords,tn,nUnitCell_1,nUnitCell_2) &
   !$omp reduction(+:zsum) 
   do i = 1 , ndim
    do j = 1 , ndim

    icount=1 
    
    rij = norm2(Coords(i,:)-Coords(j,:))

    zsum = zsum - 0.5_dp*fv(Coords(i,:),Coords(j,:))*conjg(zFock(i,j,icount))*zFock(i,j,icount)
  
    do icount=2,numNeighborCells
      n1=nUnitCell_1(icount)
      n2=nUnitCell_2(icount)

      zsum = zsum - 0.5_dp*fv(Coords(i,:),[Coords(j,1:2)-n1*tn(:,1)-n2*tn(:,2),Coords(j,3)])*conjg(zFock(i,j,icount))*zFock(i,j,icount)
      
      zsum=zsum - 0.5_dp*fv(Coords(i,:),[Coords(j,1:2)+n1*tn(:,1)+n2*tn(:,2),Coords(j,3)])*zFock(j,i,icount)*conjg(zFock(j,i,icount))
    
    enddo

    end do
  end do
  !$omp end parallel do


  fsum=real(zsum,dp)


end subroutine GetFockEnergy


subroutine OptimalStep(Step,HartreeEnergy,FockEnergy,HubbardEnergy,KineticEnergy,zFock,zFockIn,&
  Density,DensityIn,DensitySub,ndim,numNeighborCells,Coords,LongRange,Delta,nUnitCell_1,nUnitCell_2,tn)

  integer(dp) , intent(in) :: ndim, numNeighborCells
  integer(dp) , intent(in) :: nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells)
  complex(dp), intent(in) :: zFock(ndim,ndim,numNeighborCells,numS), zFockIn(ndim,ndim,numNeighborCells,numS)
  real(dp), intent(in) :: Coords(ndim,3), LongRange(ndim,ndim), tn(2,3), Delta
  real(dp), intent(in) :: Density(ndim,numS), DensityIn(ndim,numS), DensitySub(ndim)
  real(dp), intent(inout) :: HartreeEnergy,FockEnergy(numS),HubbardEnergy,KineticEnergy(numS)
  real(dp), intent(out) :: Step

  real(dp) :: denstemp(ndim,numS), temp1, temp2, coeffs2, coeffs1, coeffs0
  integer(dp) :: nspin

    coeffs0 = alpha*sum(FockEnergy)*real(3-numS,dp) + sum(KineticEnergy)*real(3-numS,dp) + U*HubbardEnergy + alphaH*HartreeEnergy

    do nspin=1,numS
      denstemp(:,nspin) = Density(:,nspin) - DensitySub(:)
    enddo

    do nspin=1,numS
      call GetKineticEnergy(KineticEnergy(nspin),Coords,zFock(:,:,:,nspin),Delta,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,tn)
      call GetFockEnergy(FockEnergy(nspin),Coords,zFock(:,:,:,nspin),nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,tn)
    enddo
    
    if(numS.EQ.1)then
      call GetHartreeHubbardEnergy(HartreeEnergy,HubbardEnergy,LongRange,denstemp(:,1),denstemp(:,1),ndim)
    else
      call GetHartreeHubbardEnergy(HartreeEnergy,HubbardEnergy,LongRange,denstemp(:,1),denstemp(:,2),ndim)
    endif

    temp1 = alpha*sum(FockEnergy)*real(3-numS,dp) + sum(KineticEnergy)*real(3-numS,dp) + U*HubbardEnergy + alphaH*HartreeEnergy

    do nspin=1,numS
      denstemp(:,nspin) = .5_dp*(Density(:,nspin) + DensityIn(:,nspin)) - DensitySub(:)
    enddo
    
    do nspin=1,numS
      call GetKineticEnergy(KineticEnergy(nspin),Coords,.5_dp*(zFockIn(:,:,:,nspin)+zFock(:,:,:,nspin)),Delta,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,tn)
      call GetFockEnergy(FockEnergy(nspin),Coords,.5_dp*(zFockIn(:,:,:,nspin)+zFock(:,:,:,nspin)),nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,tn)
    enddo
    
    if(numS.EQ.1)then
      call GetHartreeHubbardEnergy(HartreeEnergy,HubbardEnergy,LongRange,denstemp(:,1),denstemp(:,1),ndim)
    else
      call GetHartreeHubbardEnergy(HartreeEnergy,HubbardEnergy,LongRange,denstemp(:,1),denstemp(:,2),ndim)
    endif

    temp2 =  alpha*sum(FockEnergy)*real(3-numS,dp) + sum(KineticEnergy)*real(3-numS,dp) + U*HubbardEnergy + alphaH*HartreeEnergy

    coeffs1 = 4.0*temp2 - temp1 - 3.0*coeffs0
    coeffs2 = temp1 - coeffs0 - coeffs1

    ! chose 0<Step<=1 to minimize the energy (or StepAlternative if not possible)
    if((coeffs2.gt.0).and.(abs(-coeffs1/2.0_dp/coeffs2 - .5_dp).lt.0.5_dp))then
      Step = -coeffs1/2.0/coeffs2
    elseif(temp1.lt.coeffs0)then
      Step = 1.0_dp
    else
      Step = StepAlternative
    endif

    do nspin=1,numS
      denstemp(:,nspin) = (1-Step)*DensityIn(:,nspin) + Step*Density(:,nspin) - DensitySub(:)
    enddo

    do nspin=1,numS
      call GetKineticEnergy(KineticEnergy(nspin),Coords,&
        (1.0_dp-Step)*zFockIn(:,:,:,nspin)+Step*zFock(:,:,:,nspin),Delta,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,tn)
      call GetFockEnergy(FockEnergy(nspin),Coords,&
        (1.0_dp-Step)*zFockIn(:,:,:,nspin)+Step*zFock(:,:,:,nspin),nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,tn)
    enddo
    if(numS.EQ.1)then
      call GetHartreeHubbardEnergy(HartreeEnergy,HubbardEnergy,LongRange,denstemp(:,1),denstemp(:,1),ndim)
    else
      call GetHartreeHubbardEnergy(HartreeEnergy,HubbardEnergy,LongRange,denstemp(:,1),denstemp(:,2),ndim)
    endif


end subroutine OptimalStep

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure function fv(r_i,r_j)

  real(dp) , intent(in)    :: r_i(3), r_j(3)
  real(dp) :: fv, x_i, x_j, y_i, y_j, z_i, z_j, r, r2
  integer(dp) n

  x_i = r_i(1)
  y_i = r_i(2)
  z_i = r_i(3)
  x_j = r_j(1)
  y_j = r_j(2)
  z_j = r_j(3)
  r = norm2(r_i - r_j)
  r2 = sqrt((x_i-x_j)**2 + (y_i-y_j)**2)

  if(nscreen.eq.2)then
     fv=0._dp

     if(r.GT.0.01_dp)then
        if((r2.LT.xi).or.(r2.lt.2*nlayers*tz))then
          do n=-20,20
            fv  = fv + (-1.0_dp)**n/sqrt((x_j - x_i)**2 + (y_j - y_i)**2 + (z_j - (n*xi + z_i*(-1.0_dp)**n))**2)
          enddo
        else 
          fv = 2.0_dp*sqrt(2.0_dp)/sqrt(r2/xi)*exp(-pi*r2/xi)/xi
        endif

     else
        fv=0.0_dp
     endif

  endif

  if(nscreen.EQ.1)then
    if(r.GT.0.01_dp)then
      fv = 1.0_dp/r - 1.0_dp/sqrt((x_i-x_j)**2 + (y_i-y_j)**2 + (z_i + z_j + xi)**2._dp)    
    else
      fv=0.0_dp
    endif
  endif

  return
end function fv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end module HartreeFock
