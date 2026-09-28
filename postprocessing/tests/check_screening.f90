program check_screening
use Setup, only: dp, ndim, numI, numS, numC, alpha, alphaH, U, Delta, fphase, nscreen
use Geometry, only: NumberNeighborCells, OrderNeighborCells
use HartreeFock, only: LongRangeInteraction, HamiltonianHartreeFock
use Hamiltonian, only: longRange, ComplexHk_MF, vxy12I_MF
use PostProcessingInput, only: BuildGeometry, ValidateResponseModel
implicit none
integer(dp) :: nc, i, j, m, point, axis
integer(dp), allocatable :: n1(:), n2(:)
integer(dp) :: nearestUC(ndim,3), nearestT(ndim,3)
real(dp) :: coords(ndim,3),t1(2),t2(2),t3(2),g1(2),g12(2),rotation(2,2),tn(2,3)
real(dp) :: lr(ndim,ndim), lr_ref(ndim,ndim), potential(ndim),density(ndim),k(2),kp(2),km(2)
real(dp) :: h_error, v_error, hstep
complex(dp), allocatable :: rho(:,:,:)
complex(dp) :: h(ndim,ndim),ref(ndim,ndim),hp(ndim,ndim),hm(ndim,ndim)
complex(dp) :: v1(ndim,ndim),v2(ndim,ndim),vi(ndim,ndim),expected

call ValidateResponseModel()
call NumberNeighborCells(numI,nc)
allocate(n1(nc),n2(nc),rho(ndim,ndim,nc))
call OrderNeighborCells(numI,nc,n1,n2)
call BuildGeometry(coords,t1,t2,t3,g1,g12,rotation)
tn=reshape([t1,t2,t3],[2,3])
! Slater-Koster does not use nearest-neighbor tables.
nearestUC=0
nearestT=0
do m=1,nc
  do j=1,ndim
    do i=1,ndim
      rho(i,j,m)=cmplx(0.03_dp*cos(real(3*i+2*j+m,dp)), &
                         0.02_dp*sin(real(i+4*j+2*m,dp)),dp)
    enddo
  enddo
enddo
do i=1,ndim
  rho(i,i,1)=0.5_dp+0.01_dp*cos(real(i,dp))
  do j=1,i-1
    rho(j,i,1)=conjg(rho(i,j,1))
  enddo
enddo
density=real([(rho(i,i,1),i=1,ndim)],dp)-0.5_dp
call longRange(lr,coords,ndim,t1,t2)
call LongRangeInteraction(lr_ref,coords,ndim,t1,t2)
if(maxval(abs(lr-lr_ref))>1e-12_dp) error stop 'Hartree screening mismatch'
potential=2*alphaH*matmul(lr,density)+U*density
h_error=0
v_error=0
hstep=1e-6_dp
do point=1,3
  select case(point)
  case(1)
    k=0
  case(2)
    k=0.137_dp*g1+0.271_dp*g12
  case(3)
    k=0.5_dp*g1
  end select
  call ComplexHk_MF(h,coords,potential,alpha,rho,n1,n2,ndim,nc,k,tn)
  call HamiltonianHartreeFock(ref,coords,potential,alpha,Delta,rho,n1,n2,ndim,nc,k,tn,nearestUC,nearestT)
  h_error=max(h_error,maxval(abs(h-ref)))
  ! Complete the upstream lower triangle for comparison to upper-triangle currents.
  do i=1,ndim
    do j=i+1,ndim
      ref(i,j)=conjg(ref(j,i))
    enddo
  enddo
  do axis=1,2
    kp=k
    km=k
    kp(axis)=kp(axis)+hstep
    km(axis)=km(axis)-hstep
    call HamiltonianHartreeFock(hp,coords,potential,alpha,Delta,rho,n1,n2,ndim,nc,kp,tn,nearestUC,nearestT)
    call HamiltonianHartreeFock(hm,coords,potential,alpha,Delta,rho,n1,n2,ndim,nc,km,tn,nearestUC,nearestT)
    call vxy12I_MF(v1,v2,vi,coords,alpha,rho,n1,n2,ndim,nc,k,tn,axis)
    do i=1,ndim
      do j=i,ndim
        expected=-conjg(hp(j,i)-hm(j,i))/(2*hstep) &
          +cmplx(0.0_dp,real(1-fphase,dp)*(coords(i,axis)-coords(j,axis)),dp)*ref(i,j)
        v_error=max(v_error,abs(v1(i,j)+v2(i,j)+vi(i,j)-expected))
      enddo
    enddo
  enddo
enddo
write(*,'(A,I0,A,ES12.4,A,ES12.4)') 'nscreen=',nscreen,': H error=',h_error,', current error=',v_error
if(h_error>1e-10_dp) error stop 'Hamiltonian differs from dev'
if(v_error>1e-7_dp) error stop 'current differs from dev Hamiltonian derivative'
end program check_screening
