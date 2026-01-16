
! 
! gfortran -O3 -march=native -fopenmp -o bilayer LapackRoutines.f90 Hamiltonian.f90 HartreeFock_1.f90 OrderParameterPlot.f90 compute_orderparams.f90 -llapack -lblas

program HF_longrange

use omp_lib
use Hamiltonian
use HartreeFock
use OrderPlot
use lapack_routines

 implicit none

 character(250) :: filename
 character(210) :: parameters, parametersf
 character(210) :: parametersAlpha
  character(150) :: my_iomsg
  character(2) :: kOP
  character(1) :: dop
 integer(dp) :: ivk1, ivk2, ivk, n, m, pos, i,j
 integer(dp) :: flag=0_dp

 
 integer(dp) :: icount,ncount,ncount1,ncount2,nDirac,nsymBreak,ndeg1,ndeg2
 integer(dp) :: ndimEV,nNP,nvb,il,iu,nband,nspin,it,itstop,nflagMax,nflagMin,nstop
 integer(dp) :: nmu=0_dp,nmuDeg(numS),ndeg(numS)
 integer(dp) :: i1max,i1min,i2max,i2min,nBottom,nmax(2),n1,n2,iUnitCell,ind
 integer(dp) :: nChiral,nChern,nKekule,nEnergy,nConvergence
 integer(dp) :: my_iostat, rcc

 real(dp) :: ang, angulo, ayu_i, cs, sn, eDiracT, time1, time2, fmu,SB
 real(dp) :: dmaxIn(2),dmax(2),fmaxBand,fmin,sum,sumd,sumn,temp,AbsErr

 real(dp) :: vk(2), t1(2), t2(2), t3(2), vq1(2), vq12(2)
 real(dp) :: ang_mat(2,2),Gn(8,2),LatticeVectors(0:6,2), Ac=sqrt(3._dp)/2._dp*real(ndim,dp)
 integer(dp) :: FockIndexToVector(0:6,2)

 integer(dp) :: nBandsMu(2,2), nflag(2), nstate(2),maxBand(2),nBandMu,nbandMax,nbandMin
 real(dp) :: eDirac(2), dcheck(2), density(2)
 real(dp) :: Gap(2), Band(2), fmax(2)

! real(dp) :: GSEnergy(2)=0.0_dp, energy0(2)=0.0_dp, energyFS=0.0_dp
! real(dp) :: GSEnergyKinetic(numS),GSEnergyFock(numS),GSEnergyHartree,GSEnergyU
! real(dp) :: GSEnergyKineticIn(numS),GSEnergyFockIn(numS),GSEnergyHartreeIn,GSEnergyUIn

 real(dp) :: vkx,vky,det,Delta
 
 complex(dp) :: ztest,zinput,zinput1,zinput2
 real(dp) :: test

 integer(dp), allocatable :: npointsBZ(:,:), npointsBack(:,:), ink(:,:,:) !,  nsym(:,:)
 integer(dp), allocatable :: ind_j(:), ind_l(:), nn(:,:), nt(:,:)
 real(dp), allocatable    :: coord(:,:), vkBZ(:,:,:) !,fsum(:,:)
! real(dp), allocatable    :: bandsTB(:,:,:),sortEnergia(:,:),sortEnergiaBoth(:,:),  energia(:,:)
!  complex(dp), allocatable :: zh(:,:,:) !,zh0(:,:,:)
!  complex(dp), allocatable :: zevTB(:,:,:,:),zsortEV(:,:,:)

 real(dp), allocatable :: dens(:,:) !, densIn(:,:) !, dens0(:,:), densFS(:)
!  real(dp), allocatable :: potential(:,:)
!  real(dp), allocatable :: fLongRange(:,:)
 complex(dp), allocatable :: zfock(:,:,:,:),zfock0(:,:,:,:) !,zfock(:,:,:,:),zfockIn(:,:,:,:) !,zP(:,:)
 
!  real(dp) :: fmassChiralS(2),fmassChiralA(2),fmassChiralD(2),fmassChiralT(2)
!  real(dp) :: fmassPS(2),fmassPA(2),fmassSS(2),fmassSA(2)
!  real(dp) :: fmassPSAdd(2),fmassPAAdd(2),fmassSSAdd(2),fmassSAAdd(2)
!  real(dp) :: fmassPSuc,fmassPAuc,fmassSSuc,fmassSAuc
!  real(dp) :: fmassPSaltUC,fmassPAaltUC,fmassSSaltUC,fmassSAaltUC
!  real(dp) :: fmassPSFull,fmassPAFull,fmassSSFull,fmassSAFull
!  real(dp), allocatable :: densChiral(:,:),densChernSin(:,:),densChernCos(:,:),densChernAbs(:,:)
!  real(dp), allocatable :: densChernSinAdd(:,:),densChernCosAdd(:,:),densChernAbsAdd(:,:)

!  integer(dp) :: nKekuleA1,nKekuleB1,nKekuleA2,nKekuleB2,nKekuleMin,nKekuleT,nKekuleAT,nKekuleTime,nKekuleOP,nKekuleOPr,nKekuleOPi
 integer(dp), allocatable :: nnkek(:,:,:), dictkek(:,:,:)
!  integer(dp), allocatable  :: dictKekule(:),kekulePlaquete(:,:)
 real(dp), allocatable :: coordKekule(:,:)
 real(dp), allocatable :: current(:,:,:),currentFull(:,:,:)
!  complex(dp), allocatable :: zcurrentAdd(:,:,:),zcurrentOP(:,:,:),zchiral(:,:,:)
 complex(dp), allocatable :: fkbkpa(:), fkakpb(:), fkpk(:)
 real(dp), allocatable :: dens1ndim(:), dens2ndim(:)
!  real(dp) :: vc1(3),vc2(3)
!  complex(dp) :: zvct1(3),zvct2(3)
!  real(dp) :: vc1full(3),vc2full(3)
!  complex(dp) :: zvc1add(3),zvc2add(3),zvc1OP(4),zvc2OP(4)
 

 integer(dp), parameter :: LDA=ndim
 integer(dp), parameter :: LWORK=2*ndim
 integer(dp) :: INFO
 complex(dp), allocatable  :: ZWORK(:)
 real(dp), allocatable  :: RWORK(:)

 
 allocate(ZWORK(LWORK))
 allocate(RWORK(3*ndim))

 allocate(vkBZ(1:numk,1:numk,1:2))

!  allocate(densFS(1:ndim))
!  allocate(dens0(1:ndim,1:numS))
 allocate(dens(1:ndim,1:numS))
!  allocate(densIn(ndim,1:numS))
!  allocate(energia(ndim,1:numS))
 allocate(coord(ndim,3))
!  allocate(zh(ndim,ndim,1:numS))
!  allocate(zh0(1:ndim,1:ndim,1:numS))
!  allocate(zP(1:ndim,1:1))

!  allocate(coordKekule(ndim,3))
!  allocate(dictKekule(ndim))
!  allocate(kekulePlaquete(ndim,12))

 allocate(nnkek(ndim/2,2,3))
 allocate(dictKek(ndim/2,6,2))

!  allocate(current(ndim,3,numS))
!  allocate(currentFull(ndim,3,numS))
!  allocate(zcurrentAdd(ndim,3,numS))
!  allocate(zchiral(ndim,3,numS))
!  allocate(zcurrentOP(ndim,3,numS))

 allocate(nn(1:ndim,1:3))
 allocate(nt(1:ndim,1:3))

 allocate(fkbkpa(ndim/2))
 allocate(fkakpb(ndim/2))
 allocate(fkpk(ndim))
 allocate(dens1ndim(ndim))
 allocate(dens2ndim(ndim))

!  allocate(densChiral(1:ndim,1:numS))
!  allocate(densChernSin(1:ndim,1:numS))
!  allocate(densChernCos(1:ndim,1:numS))
!  allocate(densChernAbs(1:ndim,1:numS))
!  allocate(densChernSinAdd(1:ndim,1:numS))
!  allocate(densChernCosAdd(1:ndim,1:numS))
!  allocate(densChernAbsAdd(1:ndim,1:numS))

!  allocate(zsortEV(1:ndim,1:numb*numk*numk,1:numS))
!  allocate(sortEnergia(1:numb*numk*numk,1:numS))

!  allocate(nsym(6*ntheta*(ntheta+1),3))
!  if(numS.EQ.2)then
!     allocate(sortEnergiaBoth(1:2*numb*numk*numk,1:2))
!  endif

!  allocate(ink(1:numb*numk*numk,1:2,1:numS))

!  allocate(Emin(1:numb,1:numS))
!  allocate(Emax(1:numb,1:numS))
!  allocate(fsum(1:numb,1:numS))

!  allocate(potential(1:ndim,1:numS))
!  allocate(fLongRange(1:ndim,1:ndim))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  Set long-range interaction cells
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ind=0_dp

do n1=-numI,numI
do n2=-numI,numI

    if(n1+n2.LE.numI.AND.n1+n2.GE.-numI)then
    if(n1.NE.0_dp.OR.n2.NE.0_dp)then
    ind=ind+1_dp
    endif
    endif

enddo
enddo

ind=ind/2+1 ! only take half +1 due to symmetry zfock(i,j,n1,n2)=congj(zfock(j,i,-n1,-n2))


!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!

! allocate(zfockFS(1:ndim,1:ndim,1:ind))
allocate(zfock0(1:ndim,1:ndim,1:ind,1:numS))
allocate(zfock(1:ndim,1:ndim,1:ind,1:numS))
! allocate(zfockIn(1:ndim,1:ndim,1:ind,1:numS))
! zfockFS(:,:,:) =  cmplx(0.0_dp,0.0_dp,dp)
zfock0(:,:,:,:) = cmplx(0.0_dp,0.0_dp,dp)
zfock(:,:,:,:) = cmplx(0.0_dp,0.0_dp,dp)
! zfockIn(:,:,:,:) = cmplx(0.0_dp,0.0_dp,dp)

allocate(ind_j(1:ind))
allocate(ind_l(1:ind))

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!

ind_j=0_dp
ind_l=0_dp
icount=1

ind_j(icount)=0_dp
ind_l(icount)=0_dp

if(numI.GT.0)then

    icount=icount+1
    ind_j(icount)=-1_dp
    ind_l(icount)=0_dp

    icount=icount+1
    ind_j(icount)=-1_dp
    ind_l(icount)=1_dp

    icount=icount+1
    ind_j(icount)=0_dp
    ind_l(icount)=-1_dp


endif

if(numI.GT.1)then

    do n1=-numI,numI
    do n2=-numI,numI

        if(n1+n2.LE.numI.AND.n1+n2.GE.-numI.AND.icount.LT.ind)then
        if(n1.NE.0.OR.n2.NE.0)then
        icount=icount+1
        ind_j(icount)=n1
        ind_l(icount)=n2
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!  Set filename
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(nfilling.LT.0._dp)then
    dop='-'
else
    dop='+'
endif

write(parameters,'(A2,I0,A5,I0,A6,F3.1,A6,F4.2,A2,F4.2,A8,A1,F3.1,A6,I0,A6,I0,A2,A5,I0,A3,I0,A4)') '-I',ntheta,'-numk',numk,'-tperp',tperp, &
'-alpha',0.59_dp,'-U',4.0_dp,'-filling','+',0.0,'-relax',nrelax,'-SymBr',1,'IV','-numI',numI,'-dp',dp,'.dat'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Set geometry Brillouin zone
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! area of the supercell

ayu_i = 3.0_dp*ntheta**2 + 3.0_dp*ntheta + 1.0_dp

! reciprocal lattice vectors of the superlattice

vq1  = vkd/ayu_i*(real(3*ntheta+1,dp)*a1+a2)
vq12 = vkd/ayu_i*(real(3*ntheta+2,dp)*a2-a1)

! angle of rotation cs

cs = 1.0_dp-1.0_dp/(2.0_dp*ayu_i)
angulo = acos(cs)
ang    = 0.5_dp*angulo
ang_mat = reshape([cos(ang),-sin(ang),sin(ang),cos(ang)],[2,2])
sn = sqrt(1.0_dp-cs**2)

vq1  = matmul(ang_mat,vq1 )
vq12 = matmul(ang_mat,vq12)

Gn(1,:)=vq1(:)
Gn(2,:)=vq12(:)
Gn(3,:)=vq12(:)-vq1(:)
Gn(4,:)=-Gn(1,:)
Gn(5,:)=-Gn(2,:)
Gn(6,:)=-Gn(3,:)
Gn(7,:)=Gn(1,:)
Gn(8,:)=Gn(2,:)


ncount=numk*numk

allocate(npointsBZ(ncount,8))
allocate(npointsBack(numk,numk))

call samplePoints1(npointsBZ,npointsBack,numk,ncount)

do ivk2=0,numk-1
do ivk1=0,numk-1

    vkBZ(ivk1+1,ivk2+1,1)=(ivk1*vq1(1)+ivk2*vq12(1))/real(numk,dp)
    vkBZ(ivk1+1,ivk2+1,2)=(ivk1*vq1(2)+ivk2*vq12(2))/real(numk,dp)

enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Set geometry real space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Moire lattice parameters

t1 = real( ntheta  ,dp)*a1 + real(  ntheta+1,dp)*a2
t2 = real(-ntheta-1,dp)*a1 + real(2*ntheta+1,dp)*a2
t3 = t2-t1


call WignerSeitzC6(coord,t1,t2,t3,cs,sn)

!!!!!! Rotate cell to symmetrize around x-axis
do n=1,ndim
    coord(n,:) = matmul(ang_mat,coord(n,1:2))
end do


t1 = matmul(ang_mat,t1)
t2 = matmul(ang_mat,t2)
t3 = t2-t1
 
LatticeVectors(0,:)=[0.0_dp,0.0_dp]
LatticeVectors(1,:)=t1
LatticeVectors(2,:)=t2
LatticeVectors(3,:)=t3
LatticeVectors(4,:)=-t1
LatticeVectors(5,:)=-t2
LatticeVectors(6,:)=-t3

FockIndexToVector(0,1) = 0
FockIndexToVector(0,2) = 0
FockIndexToVector(1,1) = 1
FockIndexToVector(1,2) = 0
FockIndexToVector(2,1) = 0
FockIndexToVector(2,2) = 1
FockIndexToVector(3,1) = -1
FockIndexToVector(3,2) = 1
FockIndexToVector(4,1) = -1
FockIndexToVector(4,2) = 0
FockIndexToVector(5,1) = 0
FockIndexToVector(5,2) = -1
FockIndexToVector(6,1) = 1
FockIndexToVector(6,2) = -1
 

if(nrelax.EQ.1)then
     call Relaxation(coord,vq1,vq12)
endif
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Kekule environment !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call nearestNeighbour(nn,nt,coord,ndim,LatticeVectors)
call kekulelattice(coord,ndim,ang_mat,a1,a2,LatticeVectors,dictkek)
call kekulenn(nnkek,coord,ndim,a1,a2,ang_mat,LatticeVectors)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! read zfock !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do nspin=1,numS
       
    write(filename,'(A33,A10,I0,A80)') dirFock,'Fock-nspin',nspin,parameters
    open(12, file=filename , form='unformatted' , status='old',access='direct',recl=dp*2)
    rcc=0
    do i=1,ndim
       do j=1,i-1
       !do j=1,ndim
           rcc = rcc+1
           read(12,rec=rcc) zinput   
           zfock(i,j,1,nspin) = zinput
           zfock(j,i,1,nspin) = conjg(zinput)
       enddo
       rcc = rcc+1
       read(12,rec=rcc) zinput   
       zfock(i,i,1,nspin) = zinput
    enddo
     do m=2,ind
         do i=1,ndim
             do j=1,ndim
                 rcc = rcc+1
                 read(12,rec=rcc) zinput   
                 zfock(i,j,m,nspin) = zinput
             enddo
         enddo
     enddo 
    close(12)

    !    write(filename,'(A37,A11,I0,A80)') dirFockN,'Fockf-nspin',nspin,parameters
    !    open(12, file=filename , form='unformatted' , status='old',access='direct',recl=dp*2)
    !    rcc=0
    !    do i=1,ndim
    !        do j=1,i-1
    !            rcc = rcc+1
    !            read(12,rec=rcc) zinput   
    !            zfock0(i,j,1,nspin) = zinput
    !            zfock0(j,i,1,nspin) = conjg(zinput)
    !        enddo
    !        rcc = rcc+1
    !        read(12,rec=rcc) zinput   
    !        zfock0(i,i,1,nspin) = zinput
    !    enddo
    !     do m=2,ind
    !         do i=1,ndim
    !             do j=1,ndim
    !                 rcc = rcc+1
    !                 read(12,rec=rcc) zinput   
    !                 zfock0(i,j,m,nspin) = zinput
    !             enddo
    !         enddo
    !     enddo 
    !    close(12)

    ! write(filename,'(A37,A10,I0,A80)') dirFockN,'Fock-nspin',nspin,parameters
    ! open(12, file=filename , form='unformatted' , status='old',access='direct',recl=dp*2)    
    ! write(filename,'(A37,A11,I0,A80)') dirFockN,'Fockf-nspin',nspin,parameters
    ! open(13, file=filename , form='unformatted' , status='old',access='direct',recl=dp*2)
    ! rcc=0
    ! do i=1,ndim
    !     do j=1,i-1
    !         rcc = rcc+1
    !         read(12,rec=rcc) zinput1   
    !         read(13,rec=rcc) zinput2   
    !         zfock(i,j,1,nspin) = zinput1-zinput2
    !         zfock(j,i,1,nspin) = conjg(zinput1-zinput2)
    !     enddo
    !     rcc = rcc+1
    !     read(12,rec=rcc) zinput1
    !     read(13,rec=rcc) zinput2   
    !     zfock(i,i,1,nspin) = zinput1 - zinput2
    ! enddo
    ! ! do m=2,ind
    ! !     do i=1,ndim
    ! !         do j=1,ndim
    ! !             rcc = rcc+1
    ! !             read(12,rec=rcc) zinput   
    ! !             zfock(i,j,m,nspin) = zinput
    ! !         enddo
    ! !     enddo
    ! ! enddo 
    ! close(12)
    ! close(13)

 enddo


!do nspin=1,numS
!
!
!    call InterSubInterValAlt(fkakpb,fkbkpa,dictkek,nnkek,ndim,ind,ind_j,ind_l,FockIndexToVector,zfock0(:,:,:,nspin))
!    ! write(filename,'(A54,A22,I0,A80)') dir,'InterSubInterVal-nspin',nspin,parameters
!    write(filename,'(A8,A23,I0,A80)') dir,'InterSubInterValF-nspin',nspin,parameters
!    open(98,file=filename,status='replace')
!    do i=1,ndim/2
!        write(98,fmt='(E12.5,3X,E12.5,3X,E12.5,3X,E12.5)') real(fkakpb(i)), aimag(fkakpb(i)), real(fkbkpa(i)), aimag(fkbkpa(i))
!    enddo
!    close(98)
!
!    call IntraSubInterValAlt(fkpk,ndim,ind,ind_j,ind_l,nn,nt,FockIndexToVector,zfock0(:,:,:,nspin))
!    ! write(filename,'(A54,A22,I0,A80)') dir,'IntraSubInterVal-nspin',nspin,parameters
!    write(filename,'(A8,A23,I0,A80)') dir,'IntraSubInterValF-nspin',nspin,parameters
!    open(98,file=filename,status='replace')
!    do i=1,ndim
!        write(98,*) real(fkpk(i)), aimag(fkpk(i))
!    enddo
!    close(98)
!
!    call InterSubIntraValAlt(fkakpb,fkbkpa,dictkek,nnkek,ndim,ind,ind_j,ind_l,FockIndexToVector,zfock0(:,:,:,nspin))
!    ! write(filename,'(A54,A22,I0,A80)') dir,'InterSubIntraVal-nspin',nspin,parameters
!    write(filename,'(A8,A23,I0,A80)') dir,'InterSubIntraValF-nspin',nspin,parameters
!    open(98,file=filename,status='replace')
!    do i=1,ndim/2
!        write(98,*) real(fkakpb(i)), aimag(fkakpb(i))
!    enddo
!    do i=1,ndim/2
!        write(98,*) real(fkbkpa(i)), aimag(fkbkpa(i))
!    enddo
!    close(98)
!
!    call IntraSubIntraValAlt(dens1ndim,dens2ndim,ndim,nn,nt,ind,ind_j,ind_l,FockIndexToVector,zfock0(:,:,:,nspin))
!    ! write(filename,'(A54,A22,I0,A80)') dir,'IntraSubIntraVal-nspin',nspin,parameters
!    write(filename,'(A8,A23,I0,A80)') dir,'IntraSubIntraValF-nspin',nspin,parameters
!    open(98,file=filename,status='replace')
!    do i=1,ndim
!        write(98,*) dens1ndim(i), dens2ndim(i)
!    enddo
!    close(98)
!
!
!enddo

do nspin=1,numS

    call InterSubInterValAlt(fkakpb,fkbkpa,dictkek,nnkek,ndim,ind,ind_j,ind_l,FockIndexToVector,zfock(:,:,:,nspin))
    write(filename,'(A8,A22,I0,A80)') dir,'InterSubInterVal-nspin',nspin,parameters
    open(98,file=filename,status='replace')
    do i=1,ndim/2
        write(98,fmt='(E12.5,3X,E12.5,3X,E12.5,3X,E12.5)') real(fkakpb(i)), aimag(fkakpb(i)), real(fkbkpa(i)), aimag(fkbkpa(i))
    enddo
    close(98)

    call IntraSubInterValAlt(fkpk,ndim,ind,ind_j,ind_l,nn,nt,FockIndexToVector,zfock(:,:,:,nspin))
    write(filename,'(A8,A22,I0,A80)') dir,'IntraSubInterVal-nspin',nspin,parameters
    open(98,file=filename,status='replace')
    do i=1,ndim
        write(98,*) real(fkpk(i)), aimag(fkpk(i))
    enddo
    close(98)

    call InterSubIntraValAlt(fkakpb,fkbkpa,dictkek,nnkek,ndim,ind,ind_j,ind_l,FockIndexToVector,zfock(:,:,:,nspin))
    write(filename,'(A8,A22,I0,A80)') dir,'InterSubIntraVal-nspin',nspin,parameters
    open(98,file=filename,status='replace')
    do i=1,ndim/2
        write(98,*) real(fkakpb(i)), aimag(fkakpb(i))
    enddo
    do i=1,ndim/2
        write(98,*) real(fkbkpa(i)), aimag(fkbkpa(i))
    enddo
    close(98)

    call IntraSubIntraValAlt(dens1ndim,dens2ndim,ndim,nn,nt,ind,ind_j,ind_l,FockIndexToVector,zfock(:,:,:,nspin))
    write(filename,'(A8,A22,I0,A80)') dir,'IntraSubIntraVal-nspin',nspin,parameters
    open(98,file=filename,status='replace')
    do i=1,ndim
        write(98,*) dens1ndim(i), dens2ndim(i)
    enddo
    close(98)

enddo

end program
