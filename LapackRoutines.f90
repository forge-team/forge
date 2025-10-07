module lapack_routines


implicit none


interface matrix_inverse

  module procedure z_matrix_inverse, c_matrix_inverse

end interface


interface vector_overlap

  module procedure z_vector_overlap, c_vector_overlap

end interface


interface vector_mul

  module procedure z_vector_mul, c_vector_mul

end interface


interface matrix_mul

  module procedure z_matrix_mul, c_matrix_mul

end interface


interface matrix_overlap

  module procedure z_matrix_overlap, c_matrix_overlap

end interface


interface matrix_elements

  module procedure z_matrix_elements_diag, c_matrix_elements_diag, z_matrix_elements, c_matrix_elements, & 
                  z_matrix_elements_diag_real, c_matrix_elements_diag_real

end interface


interface diagonalize

 module procedure zdiagRed, cdiagRed, zdiag, cdiag

end interface


interface diagonalizeComp

 module procedure zdiagComp, cdiagComp

end interface


interface unitarize

 module procedure z_unitarize, c_unitarize

end interface


contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine z_vector_overlap(A,B,C,z_factor)
 
  complex(8), intent(in)    :: A(:,:),B(:,:)
  complex(8), intent(out) :: C 
  complex(8), intent(in), optional :: z_factor

  integer(8) :: m1
  m1 = size(A(:,1))

  C=cmplx(0.0_8,0.0_8,8)

  if(present(z_factor)) then
    call zgemm('C','N',1_8,1_8,m1,z_factor,A,m1,B,m1,cmplx(1.0_8,0.0_8,8),C,1_8) 
  else
    call zgemm('C','N',1_8,1_8,m1,cmplx(1.0_8,0.0_8,8),A,m1,B,m1,cmplx(1.0_8,0.0_8,8),C,1_8)
  endif

end subroutine z_vector_overlap

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine c_vector_overlap(A,B,C,z_factor)

  complex(4), intent(in)    :: A(:,:),B(:,:)
  complex(4), intent(out) :: C
  complex(4), intent(in), optional :: z_factor

  integer(4) :: m1
  m1 = size(A(:,1))

  C=cmplx(0.0_4,0.0_4,4)

  if(present(z_factor)) then
    call cgemm('C','N',1,1,m1,z_factor,A,m1,B,m1,cmplx(1.0_4,0.0_4,4),C,1)
  else
    call cgemm('C','N',1,1,m1,cmplx(1.0_4,0.0_4,4),A,m1,B,m1,cmplx(1.0_4,0.0_4,4),C,1)
  endif

end subroutine c_vector_overlap

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine z_matrix_overlap(A,B,C,z_factor)
 
  complex(8), intent(in)    :: A(:,:),B(:,:)
  complex(8), intent(out) :: C(:,:) 
  complex(8), intent(in), optional :: z_factor

  integer(8) :: m0,m1
  m1 = size(A(:,1))
  m0 = size(A(1,:))

  C=cmplx(0.0_8,0.0_8,8)

  if(present(z_factor)) then
    call zgemm('C','N',m0,m0,m1,z_factor,A,m1,B,m1,cmplx(1.0_8,0.0_8,8),C,m0) 
  else
    call zgemm('C','N',m0,m0,m1,cmplx(1.0_8,0.0_8,8),A,m1,B,m1,cmplx(1.0_8,0.0_8,8),C,m0)
  endif



end subroutine z_matrix_overlap

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine c_matrix_overlap(A,B,C,z_factor)

  complex(4), intent(in)    :: A(:,:),B(:,:)
  complex(4), intent(out) :: C(:,:)
  complex(4), intent(in), optional :: z_factor

  integer(4) :: m0,m1
  m1 = size(A(:,1))

  m0 = size(A(1,:))

  C=cmplx(0.0_4,0.0_4,4)

  if(present(z_factor)) then
    call cgemm('C','N',m0,m0,m1,z_factor,A,m1,B,m1,cmplx(1.0_4,0.0_4,4),C,m0)
  else
    call cgemm('C','N',m0,m0,m1,cmplx(1.0_4,0.0_4,4),A,m1,B,m1,cmplx(1.0_4,0.0_4,4),C,m0)
  endif


end subroutine c_matrix_overlap


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine z_vector_mul(A,C,z_factor)
 
  complex(8), intent(in)    :: A(:,:)
  complex(8), intent(inout) :: C(:,:)
  complex(8), intent(in), optional :: z_factor

  integer(8) :: m1
  m1 = size(A(:,1))

  if(present(z_factor)) then
    call zgemm('N','C',m1,m1,1_8,z_factor,A,m1,A,m1,cmplx(1.0_8,0.0_8,8),C,m1)
  else
    call zgemm('N','C',m1,m1,1_8,cmplx(1.0_8,0.0_8,8),A,m1,A,m1,cmplx(1.0_8,0.0_8,8),C,m1)
  endif

end subroutine z_vector_mul

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine c_vector_mul(A,C,z_factor)

 complex(4), intent(in)    :: A(:,:)
 complex(4), intent(inout) :: C(:,:)
 complex(4), intent(in), optional :: z_factor

  integer(4) :: m1
  m1 = size(A(:,1))

  if(present(z_factor)) then
    call cgemm('N','C',m1,m1,1,z_factor,A,m1,A,m1,cmplx(1.0_4,0.0_4,4),C,m1)
  else
    call cgemm('N','C',m1,m1,1,cmplx(1.0_4,0.0_4,4),A,m1,A,m1,cmplx(1.0_4,0.0_4,4),C,m1)
  endif

end subroutine c_vector_mul



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine z_matrix_mul(A,B,C)

  complex(8), intent(in)    :: A(:,:), B(:,:)
  complex(8), intent(inout) :: C(:,:)

  integer(8) :: m1,n1,k1
  m1 = size(A(:,1))
  k1 = size(A(1,:))
  n1 = size(B(1,:))
  
  call zgemm('N','N',m1,n1,k1,cmplx(1.0_8,0.0_8,8),A,m1,B,k1,cmplx(0.0_8,0.0_8,8),C,m1)

end subroutine z_matrix_mul

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine c_matrix_mul(A,B,C)

  complex(4), intent(in)    :: A(:,:), B(:,:)
  complex(4), intent(inout) :: C(:,:)
 
  integer(4) :: m1,n1,k1
  m1 = size(A(:,1))
  k1 = size(A(1,:))
  n1 = size(B(1,:))
  
  call cgemm('N','N',m1,n1,k1,cmplx(1.0_4,0.0_4,4),A,m1,B,k1,cmplx(0.0_4,0.0_4,4),C,m1)

end subroutine c_matrix_mul

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine z_matrix_elements(A,B,C,i,j,k,l)

  complex(8), intent(in)    :: A(:,:), B(:,:)
  complex(8), intent(inout) :: C(:,:)
  integer(8), intent(in), optional :: i,j,k,l

  complex(8), allocatable :: D(:,:)
  integer(8) :: m1, n1, m2

  if(present(i)) then
    m1 = size(A(:,1))
    n1 = int(j+1-i,8)
    m2 = int(l+1-k,8)
    allocate(D(m1,n1))
    call zgemm('N','N',m1,n1,m1,cmplx(1.0_8,0.0_8,8),A,m1,B(:,i:j),m1,cmplx(0.0_8,0.0_8,8),D,m1)
    call zgemm('T','N',m2,n1,m1,cmplx(1.0_8,0.0_8,8),conjg(B(:,k:l)),m1,D,m1,cmplx(0.0_8,0.0_8,8),C,m2)
  else
    m1 = size(A(:,1))
    m2 = size(B(1,:))

    allocate(D(m1,m2))
    call zgemm('N','N',m1,m2,m1,cmplx(1.0_8,0.0_8,8),A,m1,B,m1,cmplx(0.0_8,0.0_8,8),D,m1)
    call zgemm('T','N',m2,m2,m1,cmplx(1.0_8,0.0_8,8),conjg(B),m1,D,m1,cmplx(0.0_8,0.0_8,8),C,m2)
    deallocate(D)
  end if

end subroutine z_matrix_elements

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine c_matrix_elements(A,B,C,i,j,k,l)

 complex(4), intent(in)    :: A(:,:), B(:,:)
 complex(4), intent(inout) :: C(:,:)
 integer(4), intent(in), optional :: i, j, k, l

 complex(4), allocatable :: D(:,:)
 integer(4) :: m1, n1, m2

 if(present(i)) then
   m1 = int(size(A(:,1)),4)
   n1 = int(j+1-i,4)
   m2 = int(l+1-k,4)
   allocate(D(m1,n1))
   call zgemm('N','N',m1,n1,m1,cmplx(1.0_4,0.0_4,4),A,m1,B(:,i:j),m1,cmplx(0.0_4,0.0_4,4),D,m1)
   call zgemm('T','N',m2,n1,m1,cmplx(1.0_4,0.0_4,4),conjg(B(:,k:l)),m1,D,m1,cmplx(0.0_4,0.0_4,4),C,m2)
 else
   m1 = size(A(:,1))
   m2 = size(B(1,:))
   allocate(D(m1,m1))
   call cgemm('N','N',m1,m2,m1,cmplx(1.0_4,0.0_4,4),A,m1,B,m1,cmplx(0.0_4,0.0_4,4),D,m1)
   call cgemm('T','N',m2,m2,m1,cmplx(1.0_4,0.0_4,4),conjg(B),m1,D,m1,cmplx(0.0_4,0.0_4,4),C,m2)
   deallocate(D)
 end if

end subroutine c_matrix_elements

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine z_matrix_elements_diag(A,B,C)

 complex(8), intent(in)    :: A(:,:), B(:,:)
 complex(8), intent(inout) :: C(:)

 complex(8), allocatable :: D(:,:)
 complex(8) :: zdotc
 integer :: m, i

 m = size(B(:,1))
 allocate(D(m,m))
 call zhemm('L','U',m,m,cmplx(1.0_8,0.0_8,8),A,m,B,m,cmplx(0.0_8,0.0_8,8),D,m)
 C = [(zdotc(m,B(:,i),1,D(:,i),1),i=1,m)]
 deallocate(D)

end subroutine z_matrix_elements_diag

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine c_matrix_elements_diag(A,B,C)

 complex(4), intent(in)    :: A(:,:), B(:,:)
 complex(4), intent(inout) :: C(:)

 complex(4), allocatable :: D(:,:) 
 complex(4) :: cdotc
 integer :: m, i

 m = size(B(:,1))
 allocate(D(m,m))
 call chemm('L','U',m,m,cmplx(1.0_4,0.0_4,4),A,m,B,m,cmplx(0.0_4,0.0_4,4),D,m)
 C = [(cdotc(m,B(:,i),1,D(:,i),1),i=1,m)]
 deallocate(D)

end subroutine c_matrix_elements_diag

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine c_matrix_elements_diag_real(A,B,C)

 complex(4), intent(in) :: A(:,:), B(:,:)
 real(4), intent(inout) :: C(:)

 complex(4), allocatable :: D(:,:) 
 complex(4) :: cdotc
 integer :: m, i

 m = size(B(:,1))
 allocate(D(m,m))
 call chemm('L','U',m,m,cmplx(1.0_4,0.0_4,4),A,m,B,m,cmplx(0.0_4,0.0_4,4),D,m)
 C = [(real(cdotc(m,B(:,i),1,D(:,i),1),4),i=1,m)]
 deallocate(D)

end subroutine c_matrix_elements_diag_real

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine z_matrix_elements_diag_real(A,B,C)

 complex(8), intent(in) :: A(:,:), B(:,:)
 real(8), intent(inout) :: C(:)

 complex(8), allocatable :: D(:,:)
 complex(8) :: zdotc
 integer :: m, i

 m = size(B(:,1))
 allocate(D(m,m))
 call zhemm('L','U',m,m,cmplx(1.0_8,0.0_8,8),A,m,B,m,cmplx(0.0_8,0.0_8,8),D,m)
 C = [(real(zdotc(m,B(:,i),1,D(:,i),1),8),i=1,m)]
 deallocate(D)

end subroutine z_matrix_elements_diag_real


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Start Diagonalization for complex(8)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine zdiagRed(A,eigenvalues,ev,nLower,nUpper)

  implicit none

  complex(8), intent(inout) :: A(:,:)
  real(8)   , intent(inout) :: eigenvalues(:)
  character(1) , intent(in) :: ev
  integer(8) , intent(in)      :: nLower , nUpper
 
  integer(4) :: info, m, lwork, lrwork, liwork,il4,iu4
  real(8)    :: vl, vu, abstol, dlamch
  integer(4), allocatable :: iwork(:), isuppz(:)
  complex(8), allocatable :: work(:),z(:,:)
  real(8)   , allocatable :: rwork(:)

  if( (ev /= 'V') .and. (ev /= 'N') ) then
   write(*,*) 'ev has to be either "V" or "N" '
   stop
  end if

  il4=nLower

  iu4=nUpper

  allocate( isuppz(2*(nUpper-nLower+1)) )
  allocate( z(size(A,1),nUpper-nLower+1) )
  allocate(  work(1) )
  allocate( rwork(1) )
  allocate( iwork(1) )

  lwork = -1
  lrwork= -1
  liwork= -1
  call zheevr(ev,'I','L',size(A,1),A,size(A,1),vl,vu,il4,iu4,dlamch('S'),m,eigenvalues,&
              z,size(A,1),isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)

  lwork  = int( work(1))
  lrwork = int(rwork(1))
  liwork = int(iwork(1))
  deallocate(work,rwork,iwork)
  allocate(work(lwork))
  allocate(rwork(lrwork))
  allocate(iwork(liwork))

  call zheevr(ev,'I','L',size(A,1),A,size(A,1),vl,vu,il4,iu4,dlamch('S'),m,eigenvalues,&
              z,size(A,1),isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)

  A(:,1:nUpper-nLower+1) = z

  if(info.ne.0) write(*,*) 'zheev failed' , info

 
end subroutine zdiagRed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Start Diagonalization for complex(4)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine cdiagRed(A,eigenvalues,ev,nLower,nUpper)

  implicit none

  complex, intent(inout) :: A(:,:)
  real   , intent(inout) :: eigenvalues(:)
  character(1) , intent(in) :: ev
  integer , intent(in)      :: nLower , nUpper
 
  integer :: info, m, lwork, lrwork, liwork
  real    :: vl, vu, abstol, dlamch
  integer, allocatable :: iwork(:), isuppz(:)
  complex, allocatable :: work(:),z(:,:)
  real   , allocatable :: rwork(:)

  if( (ev /= 'V') .and. (ev /= 'N') ) then
   write(*,*) 'ev has to be either "V" or "N" '
   stop
  end if

  allocate( isuppz(2*(nUpper-nLower+1)) )
  allocate( z(size(A,1),nUpper-nLower+1) )
  allocate(  work(1) )
  allocate( rwork(1) )
  allocate( iwork(1) )

  lwork = -1
  lrwork= -1
  liwork= -1
  call cheevr(ev,'I','L',size(A,1),A,size(A,1),vl,vu,nLower,nUpper,dlamch('S'),m,eigenvalues,&
              z,size(A,1),isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)

  lwork  = int( work(1))
  lrwork = int(rwork(1))
  liwork = int(iwork(1))
  deallocate(work,rwork,iwork)
  allocate(work(lwork))
  allocate(rwork(lrwork))
  allocate(iwork(liwork))

  call cheevr(ev,'I','L',size(A,1),A,size(A,1),vl,vu,nLower,nUpper,dlamch('S'),m,eigenvalues,&
              z,size(A,1),isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)

  A(:,1:nUpper-nLower+1) = z

  if(info.ne.0) write(*,*) 'cheev failed' , info

 
end subroutine cdiagRed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% End Diagonalization for complex(8)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine zdiagFull(A,eigenvalues,ev)

  implicit none

  complex(8), intent(inout) :: A(:,:)
  real(8)   , intent(inout) :: eigenvalues(:)
  character(1) , intent(in) :: ev
 
  integer(4) :: info, lwork
  complex(8), allocatable :: work(:)
  real(8)   , allocatable :: rwork(:)


  if( (ev /= 'V') .and. (ev /= 'N') ) then
   write(*,*) 'ev has to be either "V" or "N" '
   stop
  end if


  allocate( work(1) )
  allocate( rwork(3*size(A,1)))


  lwork = -1

  call zheev(ev,'L',size(A,1),A,size(A,1),eigenvalues,work,lwork,rwork,info)

  lwork  = int( work(1))
  deallocate(work)
  allocate(work(lwork))

  call zheev(ev,'L',size(A,1),A,size(A,1),eigenvalues,work,lwork,rwork,info)
  if(info.ne.0) write(*,*) 'zheev failed' , info

 
end subroutine zdiagFull

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Start Full Diagonalization for complex(4)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine cdiagFull(A,eigenvalues,ev)

  implicit none

  complex, intent(inout) :: A(:,:)
  real   , intent(inout) :: eigenvalues(:)
  character , intent(in) :: ev
 
  integer :: ndim,info, lwork
  complex, allocatable :: work(:)
  real   , allocatable :: rwork(:)

  if( (ev /= 'V') .and. (ev /= 'N') ) then
   write(*,*) 'ev has to be either "V" or "N" '
   stop
  end if

  allocate( work(1) )
  allocate( rwork(3*size(A,1)))

  lwork = -1

  call cheev(ev,'L',size(A,1),A,size(A,1),eigenvalues,work,lwork,rwork,info)

  lwork  = int( work(1))
  deallocate(work)
  allocate(work(lwork))


  call cheev(ev,'L',size(A,1),A,size(A,1),eigenvalues,work,lwork,rwork,info)
  if(info.ne.0) write(*,*) 'cheev failed' , info

 
end subroutine cdiagFull

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% End Full Diagonalization for complex(4)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Diagonalization for complex(4)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine cdiag(A,eigenvalues,ev,d)

  implicit none

  complex(4), intent(inout) :: A(:,:)
  real(4)   , intent(inout):: eigenvalues(:)
  character(1) , intent(in) :: ev
  character(1) , optional , intent(in) :: d
  
  integer(4) :: info
  integer(4), allocatable :: iwork(:)
  complex(4), allocatable :: work(:)
  real(4)   , allocatable :: rwork(:)

 
  if( (ev /= 'V') .and. (ev /= 'N') ) then
   write(*,*) 'ev has to be either "V" or "N" '
   stop
  end if

  if(present(d)) then

    if(ev=='V') then
      allocate( rwork(1+5*size(A,1)+2*size(A,1)**2) )
      allocate(  work(2*size(A,1)+size(A,1)**2) )
      allocate( iwork(3+5*size(A,1)) )
    else
      allocate( rwork(size(A,1)) )
      allocate(  work(size(A,1)+1) )
      allocate( iwork(1) )
    end if

    call cheevd(ev, 'L', size(A,1), A, size(A,1), eigenvalues, work, size(work), rwork, size(rwork), iwork ,size(iwork), info)

    if(info.ne.0) write(*,*) 'cheevd failed' , info

  else

    allocate( rwork(3*size(A,1)) )
    allocate(  work(2*size(A,1)) )

    call cheev(ev, 'L', size(A,1), A, size(A,1), eigenvalues, work, 2*size(A,1), rwork, info)

    if(info.ne.0) write(*,*) 'cheev failed' , info

  end if
  
end subroutine cdiag


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Diagonalization for complex(8)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine zdiag(A,eigenvalues,ev,d)

  implicit none

  complex(8), intent(inout) :: A(:,:)
  real(8)   , intent(inout) :: eigenvalues(:)
  character(1) , intent(in) :: ev
  character(1) , optional , intent(in) :: d
  
  integer(4) :: info
  integer(4), allocatable :: iwork(:)
  complex(8), allocatable :: work(:)
  real(8)   , allocatable :: rwork(:)

  if( (ev /= 'V') .and. (ev /= 'N') ) then
   write(*,*) 'ev has to be either "V" or "N" '
   stop
  end if

  if(present(d)) then

    if(ev=='V') then
      allocate( rwork(1+5*size(A,1)+2*size(A,1)**2) )
      allocate(  work(2*size(A,1)+size(A,1)**2) )
      allocate( iwork(3+5*size(A,1)) )
    else
      allocate( rwork(size(A,1)) )
      allocate(  work(size(A,1)+1) )
      allocate( iwork(1) )
    end if

       call zheevd(ev, 'L', size(A,1), A, size(A,1), eigenvalues, &
                   work, size(work), rwork, size(rwork), iwork ,size(iwork), info)

    if(info.ne.0) write(*,*) 'zheevd failed' , info

  else

    allocate( rwork(3*size(A,1)) )
    allocate(  work(2*size(A,1)) )

    call zheev(ev, 'L', size(A,1), A, size(A,1), eigenvalues, work, 2*size(A,1), rwork, info)

    if(info.ne.0) write(*,*) 'zheev failed' , info

  end if
  
end subroutine zdiag


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% End Diagonalization for complex(8)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Diagonalization of general complex matrix
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine zdiagComp(A,eigenvalues,ev)

  implicit none

  complex(8), intent(inout) :: A(:,:)
  complex(8)   , intent(inout) :: eigenvalues(:)
  character(1) , intent(in) :: ev
 
  integer(4) :: info, lwork
  complex(8), allocatable :: work(:)
  real(8)   , allocatable :: rwork(:)

  complex(8) , allocatable :: vl(:,:),vr(:,:)

  if( (ev /= 'V') .and. (ev /= 'N') ) then
   write(*,*) 'ev has to be either "V" or "N" '
   stop
  end if

  allocate( work(1) )
  allocate( rwork(3*size(A,1)))

  allocate( vl(size(A,1),size(A,1)))
  allocate( vr(size(A,1),size(A,1)))

  lwork = -1

  call zgeev(ev,ev,size(A,1),A,size(A,1),eigenvalues,vl,size(A,1),vr,size(A,1),work,lwork,rwork,info)

  lwork  = int( work(1))
  deallocate(work)
  allocate(work(lwork))

     call zgeev(ev,ev,size(A,1),A,size(A,1),eigenvalues,vl,size(A,1),vr,size(A,1),work,lwork,rwork,info)
  if(info.ne.0) write(*,*) 'zheev failed' , info

 
end subroutine zdiagComp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Start Full Complex Diagonalization for complex(4)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine cdiagComp(A,eigenvalues,ev)

  implicit none

  complex(4), intent(inout) :: A(:,:)
  complex(4)   , intent(inout) :: eigenvalues(:)
  character(1) , intent(in) :: ev
 
  integer(4) :: info, lwork
  complex(4), allocatable :: work(:)
  real(4)   , allocatable :: rwork(:)

  complex(4) , allocatable :: vl(:,:),vr(:,:)

  if( (ev /= 'V') .and. (ev /= 'N') ) then
   write(*,*) 'ev has to be either "V" or "N" '
   stop
  end if

  allocate( work(1) )
  allocate( rwork(3*size(A,1)))

  allocate( vl(size(A,1),size(A,1)))
  allocate( vr(size(A,1),size(A,1)))

  lwork = -1

  call cgeev(ev,ev,size(A,1),A,size(A,1),eigenvalues,vl,size(A,1),vr,size(A,1),work,lwork,rwork,info)

  lwork  = int( work(1))
  deallocate(work)
  allocate(work(lwork))


  call cgeev(ev,ev,size(A,1),A,size(A,1),eigenvalues,vl,size(A,1),vr,size(A,1),work,lwork,rwork,info)
  if(info.ne.0) write(*,*) 'cheev failed' , info

 
end subroutine cdiagComp


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Inverse of complex matrix double precision
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine z_matrix_inverse(A)

  implicit none

  complex(8), intent(inout) :: A(:,:)
 
  integer(4) :: info,n,m
  integer(8) , allocatable :: ipiv(:)
  complex(8) , allocatable :: work(:)


  allocate(ipiv(size(A,1)))
  allocate(work(size(A,1)*size(A,1)))

  call ZGETRF(size(A,1), size(A,1), A, size(A,1), ipiv, info)
  
  if(info.ne.0) write(*,*) 'zgetrf failed' , info

  call ZGETRI(size(A,1), A, size(A,1), ipiv, work, size(A,1)*size(A,1), info)

  if(info.ne.0) write(*,*) 'zgetri failed' , info


end subroutine z_matrix_inverse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Inverse of complex matrix single precision
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine c_matrix_inverse(A)

  implicit none

  complex(4), intent(inout) :: A(:,:)
 
  integer(4) :: info,n,m
  integer(4) , allocatable :: ipiv(:)
  complex(4) , allocatable :: work(:)


  allocate(ipiv(size(A,1)))
  allocate(work(size(A,1)*size(A,1)))

  call CGETRF(size(A,1), size(A,1), A, size(A,1), ipiv, info)
  
  if(info.ne.0) write(*,*) 'zgetrf failed' , info

  call CGETRI(size(A,1), A, size(A,1), ipiv, work, size(A,1)*size(A,1), info)

  if(info.ne.0) write(*,*) 'zgetri failed' , info

end subroutine c_matrix_inverse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Unitarize Hermitian complex matrix double precision
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine z_unitarize(A)
 
  complex(8), intent(inout) :: A(:,:)
 
  real(8), allocatable :: eigenvalues(:)
  complex(8), allocatable :: B(:,:), C(:,:)
  integer(8) :: m1, n
  m1 = size(A(:,1))

  allocate(eigenvalues(m1))
  allocate(B(m1,m1))
  allocate(C(m1,m1))

  call zdiagFull(A,eigenvalues,'V')

  do n=1,m1
      eigenvalues(n) = sign(1.0_8,eigenvalues(n))
  enddo

  B(:,:) = cmplx(0.0_8,0.0_8)
  do n=1,m1
    B(n,n) = cmplx(eigenvalues(n),0.0_8)
  enddo

  A = conjg(transpose(A))
  call z_matrix_elements(B,A,C)
  A(:,:) = C(:,:)

end subroutine z_unitarize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Unitarize Hermitian complex matrix single precision
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine c_unitarize(A)
 
  complex(4), intent(inout) :: A(:,:)
 
  real(4), allocatable :: eigenvalues(:)
  complex(4), allocatable :: B(:,:), C(:,:)
  integer(4) :: m1, n
  m1 = size(A(:,1))

  allocate(eigenvalues(m1))
  allocate(B(m1,m1))
  allocate(C(m1,m1))

  call cdiagFull(A,eigenvalues,'V')

  do n=1,m1
    eigenvalues(n) = sign(1.0_4,eigenvalues(n))
  enddo

  B(:,:) = cmplx(0.0_4,0.0_4)
  do n=1,m1
    B(n,n) = cmplx(eigenvalues(n),0.0_4)
  enddo

  A = conjg(transpose(A))
  call c_matrix_elements(B,A,C)
  A(:,:) = C(:,:)

end subroutine c_unitarize



end module
