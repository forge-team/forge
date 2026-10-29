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

subroutine HamiltonianHartreeFock(zH,Coords,Potential,alpha,Delta,zFock,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,vk,tn,NearestNeighborsUC,NearestNeighborsT)

  integer(dp), intent(in)    :: ndim, numNeighborCells

  integer(dp), intent(in)    :: nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells),NearestNeighborsUC(ndim,3),NearestNeighborsT(ndim,3)
  real(dp)   , intent(in)    :: tn(2,3), vk(2), Coords(ndim,3)
  real(dp)   , intent(in)    :: alpha, Delta, Potential(ndim)
  complex(dp) , intent(in)   :: zFock(ndim,ndim,numNeighborCells)
  complex(dp), intent(out) :: zH(ndim,ndim)

  integer(dp) :: i, j, icount, n1, n2
  real(dp)    :: rij, Coords_diff, tm(0:6,2)
  complex(dp) :: zi_vk_tn,zsum,zphase

  tm(0,:) = [0.0_dp, 0.0_dp]
  tm(1,:) = tn(:,1)
  tm(2,:) = tn(:,2)
  tm(3,:) = tn(:,3)
  tm(4,:) = -tn(:,1)
  tm(5,:) = -tn(:,2)
  tm(6,:) = -tn(:,3)

  zH(:,:) = cmplx(0.0_dp,0.0_dp,dp)

  do i = 1 , ndim
    ! Hartree part
    zH(i,i) = zH(i,i) + Potential(i)
    do j = 1 , i

      ! Kinetic term
      zphase=exp(cmplx(0.0_dp,-dot_product(vk,Coords(i,1:2)-Coords(j,1:2))*real(fphase,dp),dp))

      zH(i,j)=zH(i,j)-alpha*fv(Coords(i,:),Coords(j,:))*conjg(zFock(i,j,1))*zphase
     
      zH(i,j) = zH(i,j)+fTB(TBFunction,Coords(i,:),Coords(j,:),i,j,NearestNeighborsUC,NearestNeighborsT,Coords,tm)*zphase

      do icount=2,numNeighborCells
        n1=nUnitCell_1(icount)
        n2=nUnitCell_2(icount)
              
        zi_vk_tn = cmplx(0.0_dp,dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
        zH(i,j)=zH(i,j)+fTB(TBFunction,Coords(i,:), [Coords(j,1:2)-n1*tn(:,1)-n2*tn(:,2),Coords(j,3)],i,j,NearestNeighborsUC,NearestNeighborsT,Coords,tm)*exp(-zi_vk_tn)*zphase
        
        zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,n1*tn(:,1)+n2*tn(:,2)),dp)
        zH(i,j)=zH(i,j)+fTB(TBFunction,Coords(i,:), [Coords(j,1:2)+n1*tn(:,1)+n2*tn(:,2),Coords(j,3)],i,j,NearestNeighborsUC,NearestNeighborsT,Coords,tm)*exp(-zi_vk_tn)*zphase

  
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
  Coords,MomentaValues,nMomentaComponents,nLower,nUpper,ndim,RotateLayers,numk,Nk,numNeighborCells,nUnitCell_1,nUnitCell_2,tn,g1,g12,nC3pairs,NearestNeighborsUC,NearestNeighborsT)

  integer(dp), intent(in)    :: ndim, numk, Nk, nLower, nUpper,numNeighborCells, nC3pairs(ndim), RotateLayers(nlayers)
  integer(dp), intent(in)    :: nMomentaComponents(Nk,2),NearestNeighborsUC(ndim,3),NearestNeighborsT(ndim,3)
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
  !$omp shared(Coords,Nk,MomentaMod_C3,tn,nLower,nUpper,Bands,zEigenvectors,nUnitCell_1,nUnitCell_2,numNeighborCells,nC3pairs,g1,g12,NearestNeighborsUC,NearestNeighborsT) &
  !$omp shared(nMomentaComponents,MomentaValues,zFock) &
  !$omp reduction(+:zFockBulk)
  do icount=1,NkMod_C3 ! loop over the NkMod_C3 momenta

    icountk = MomentaMod_C3(icount)
    ivk1=nMomentaComponents(icountk,1)
    ivk2=nMomentaComponents(icountk,2)
    vk(:)=MomentaValues(ivk1+1,ivk2+1,:)

    write(*,*) icount, icountk

    zFockBulkTemp=cmplx(0.0_dp,0.0_dp,dp)

    call HamiltonianHartreeFock(zH_FockBulk1Temp,Coords,Potential,alpha,Delta,zFock,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,vk,tn,NearestNeighborsUC,NearestNeighborsT)

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
  Coords,MomentaValues,nMomentaComponents,nLower,nUpper,ndim,numk,Nk,numNeighborCells,nUnitCell_1,nUnitCell_2,tn,NearestNeighborsUC,NearestNeighborsT)

  integer(dp), intent(in)    :: ndim, numk, Nk, nLower, nUpper,numNeighborCells
  integer(dp), intent(in)    :: nMomentaComponents(Nk,2),NearestNeighborsUC(ndim,3),NearestNeighborsT(ndim,3)
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
  !$omp shared(Coords,Nk,tn,nLower,nUpper,Bands,zEigenvectors,nUnitCell_1,nUnitCell_2,numNeighborCells,NearestNeighborsUC,NearestNeighborsT) &
  !$omp shared(nMomentaComponents,MomentaValues,zFock) &
  !$omp reduction(+:zFockBulk)
  do icount=1,Nk

      write(*,*) icount

      zFockBulkTemp=cmplx(0.0_dp,0.0_dp,dp)

      ivk1=nMomentaComponents(icount,1)
      ivk2=nMomentaComponents(icount,2)
      vk(:)=MomentaValues(ivk1+1,ivk2+1,:)

      call HamiltonianHartreeFock(H,Coords,Potential,alpha,Delta,zFock,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,vk,tn,NearestNeighborsUC,NearestNeighborsT)

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
  zdeg=cmplx(DegFactor,0.0_dp,dp)

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
  zdeg=cmplx(DegFactor,0.0_dp,dp)
  
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

subroutine GetKineticEnergy(fsum,Coords,zFock,Delta,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,tn,NearestNeighborsUC,NearestNeighborsT)

  integer(dp), intent(in)    :: ndim, numNeighborCells

  integer(dp), intent(in)    :: nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells),NearestNeighborsUC(ndim,3),NearestNeighborsT(ndim,3)
  real(dp)   , intent(in)    :: tn(2,3), Coords(ndim,3), Delta
  complex(dp) , intent(in)   :: zFock(ndim,ndim,numNeighborCells)
  real(dp), intent(out) :: fsum

  integer(dp) :: i, j, icount, n1, n2
  complex(dp) :: zsum
  real(dp) :: tm(0:6,2)

  tm(0,:) = [0.0_dp, 0.0_dp]
  tm(1,:) = tn(:,1)
  tm(2,:) = tn(:,2)
  tm(3,:) = tn(:,3)
  tm(4,:) = -tn(:,1)
  tm(5,:) = -tn(:,2)
  tm(6,:) = -tn(:,3)

  zsum=cmplx(0.0_dp,0.0_dp,dp)

  !$omp parallel do &
  !$omp private(i,j,icount,n1,n2) &
  !$omp shared(ndim,zFock,numNeighborCells,Coords,tn,nUnitCell_1,nUnitCell_2,NearestNeighborsUC,NearestNeighborsT) &
  !$omp reduction(+:zsum) 
  do i = 1,ndim
   do j = 1,ndim
   
     icount=1
      
     zsum=zsum+fTB(TBFunction,Coords(i,:),Coords(j,:),i,j,NearestNeighborsUC,NearestNeighborsT,Coords,tm)*zFock(i,j,1)

      do icount=2,numNeighborCells
        n1=nUnitCell_1(icount)
        n2=nUnitCell_2(icount)
        
      zsum=zsum+fTB(TBFunction,Coords(i,:), [Coords(j,1:2)-n1*tn(:,1)-n2*tn(:,2),Coords(j,3)],i,j,NearestNeighborsUC,NearestNeighborsT,Coords,tm)*zFock(i,j,icount)
      zsum=zsum+fTB(TBFunction,Coords(i,:), [Coords(j,1:2)+n1*tn(:,1)+n2*tn(:,2),Coords(j,3)],i,j,NearestNeighborsUC,NearestNeighborsT,Coords,tm)*conjg(zFock(j,i,icount))

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
  Density,DensityIn,DensitySub,ndim,numNeighborCells,Coords,LongRange,Delta,nUnitCell_1,nUnitCell_2,tn,NearestNeighborsUC,NearestNeighborsT)

  integer(dp) , intent(in) :: ndim, numNeighborCells,NearestNeighborsUC(ndim,3),NearestNeighborsT(ndim,3)
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
      call GetKineticEnergy(KineticEnergy(nspin),Coords,zFock(:,:,:,nspin),Delta,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,tn,NearestNeighborsUC,NearestNeighborsT)
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
      call GetKineticEnergy(KineticEnergy(nspin),Coords,.5_dp*(zFockIn(:,:,:,nspin)+zFock(:,:,:,nspin)),Delta,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,tn,NearestNeighborsUC,NearestNeighborsT)
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
        (1.0_dp-Step)*zFockIn(:,:,:,nspin)+Step*zFock(:,:,:,nspin),Delta,nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,tn,NearestNeighborsUC,NearestNeighborsT)
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
        if((r2.LT.xi/2.0_dp).or.(r2.lt.2*nlayers*tz))then
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


! Postprocessing compatibility namespace (shares the dev physical settings).
! Compatibility routines grouped by the same responsibilities as FORGE dev.
module ResponseHartreeFock
use ResponseSetup
use lapack_routines
use HartreeFock, only: fv, LongRangeInteraction
use ResponseTightBinding
implicit none
contains


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
end module ResponseHartreeFock
