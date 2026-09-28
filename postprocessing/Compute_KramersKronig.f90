! Cleaned from plotKramersKronig.f90; original modules retained.
! Reads the conductivity output and applies the Kramers-Kronig transform.
! Geometry/input: FORGE dev; response kernels: Hamiltonian.f90. See README.md.

program compute_kramers_kronig

use omp_lib
use Hamiltonian
use PostProcessingInput, only: BuildGeometry, ReadFock, InputParameters, ValidateResponseModel
use Geometry, only: NumberNeighborCells, OrderNeighborCells
use Setup, only: alphaInput => alpha, alphaHInput => alphaH
use lapack_routines

 implicit none

 character(1024) :: filename
 character(220) :: parameters
 integer(dp) :: n, imu

 integer(dp) :: numkPlot, nKK, nstep, nstepI
 integer(dp) :: it

 real(dp) :: fmu
 real(dp) :: CBand
 real(dp) :: DiamagX(nmu,nT), DiamagY(nmu,nT)

 real(dp) :: om

 real(dp) :: temp1, temp2

 integer(dp) :: nspacing
 integer(dp) :: nomega, nshift
 real(dp), parameter :: Ac=4._dp*pi*pi ! normalization

 real(dp), allocatable :: sigmaTotX(:,:,:)
 real(dp), allocatable :: ReSigmaX(:,:,:)

 real(dp) :: alpha

write(*,*) 'Test',dir

! Consume the spectral filenames generated from the same FORGE state.
call InputParameters(parameters)

! Start Read

write(filename,'(A,A5,A)') trim(dir),'CBand',trim(parameters)
open(87,file=filename,status='old', action='read')

read(87,*) CBand,numkPlot,nspacing

close(87)

! numkPlot was read from the spectral output header.

do imu=1,nmu

fmu=CBand*(nmu-imu)/real(nmu,8)

write(*,*) imu,fmu

do iT=1,nT ! iT

write(filename,'(A,A5,I0,A3,I0,A9,I0,A9,I0,A)') trim(dir),'DiaXY',int(ftemperature(iT)),'-Mu',imu,'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(67,file=filename,status='old', action='read')

read(67,*) nomega,nshift,temp1,temp2

DiamagX(imu,iT)=temp1

DiamagY(imu,iT)=temp2

enddo

enddo

  allocate(sigmaTotX(nomega,nmu,nT))
  allocate(ReSigmaX(nomega,nmu,nT))

do imu=1,nmu

fmu=CBand*(nmu-imu)/real(nmu,8)

do iT=1,nT ! iT

write(*,*) imu,DiamagX(imu,iT)
write(*,*) imu,DiamagY(imu,iT)

write(filename,'(A,A15,I0,A3,I0,A9,I0,A9,I0,A)') trim(dir),'ConductivityTot',int(ftemperature(iT)),'-Mu',imu,'-numkPlot',numkPlot,'-nspacing',nspacing,trim(parameters)
open(97,file=filename,status='old', action='read')

! Normalization to "sigma_0=e^2/hbar/8"
! Implementation with "exact" symmetry for negative frequencies

do n=nshift,nomega

read(97,*) om,temp1

sigmaTotX(n,imu,iT)=temp1/(8_dp/Ac)

enddo
close(97)

enddo !iT

enddo !imu

write(*,*) 'Finish Read'

nstep=1_dp

nstepI=1_dp

nKK=min(int(nomega*0.05)+nshift,nomega)

write(*,*) 'Write nKK'

write(*,*) nKK,nshift

do imu=1,nmu

do iT=1,nT ! iT

call KramersKronig(ReSigmaX(:,imu,iT),sigmaTotX(:,imu,iT),nKK,nomega,nspacing,nshift,nstep,nstepI)

write(*,*) '1',maxval(abs(ReSigmaX(:,imu,iT)))

enddo ! iT

enddo ! imu

write(*,*) 'Output'

write(*,*) 'Test',numkPlot

do imu=1,nmu

fmu=CBand*(nmu-imu)/real(nmu,8)

write(*,*) imu,fmu

do iT=1,nT ! iT

write(filename,'(A,A7,I0,A3,I0,A9,I0,A9,I0,A6,I0,A7,I0,A4,I0,A)') trim(dir),'RealTot',int(ftemperature(iT)),'-Mu',imu,'-numkPlot',numkPlot,'-nspacing',nspacing,'-nstep',nstep,'-nstepI',nstepI,'-nKK',nKK,trim(parameters)
open(87,file=filename,status='replace')

! Normalization to "sigma_0=e^2/hbar/8"
! Implementation with "exact" symmetry for negative frequencies
do n=nshift,nKK,nstep
om=real(n-nshift)/real(nspacing)
write(87,*) om,ReSigmaX(n,imu,iT)/Ac+pi/2._dp*(DiamagX(imu,iT)+DiamagY(imu,iT))/2._dp!,(DiamagX(imu,iT)+DiamagY(imu,iT))/2._dp
enddo
close(87)

enddo !iT

enddo !imu

write(*,*) 'TestEnd',numkPlot

end program compute_kramers_kronig
