
module Setup

implicit none

! "model" parameters
integer, parameter :: dp = 8                                                               ! dp=4/8: single/double precision
integer(dp), parameter :: nlayers = 2                                                      ! number of layers
integer(dp), parameter :: ntheta = 9                                                      ! twist angle = acos(1 - 1/(6ntheta^2 + 6ntheta + 2))
integer(dp), parameter :: RotateLayers(nlayers) = [-1,+1]                                  ! Positive/negative rotation of layers 
integer(dp), parameter :: ndim = nlayers*2*(ntheta**2+(ntheta+1)**2+(ntheta+1)*ntheta)     ! number of atoms in unt cell
integer(dp), parameter :: numk = 6                      ! number of grid points in the BZ  Nk = numk * numk
integer(dp), parameter :: numkPostProc = 60             ! number of grid points in the BZ for post-processing Nk = numkPostProc * numkPostProc
integer(dp), parameter :: TBFunction = 1                 ! 1/2: Slater-Koster/Wannier tight-binding function (Wannier only for bilayer)
integer(dp), parameter :: nrelax = 0                     ! 0/1/2: unrelaxed/relaxed lattice Nam etal/relaxed lattice Carr et al (1/2 only for bilayer) 

! "solution" parameters
real(dp), parameter :: nfilling = 0.0_dp                 ! filling in units of e/unit cell. -4/+4:empty/full flat bands, 0:charge neutrality
integer(dp), parameter :: numS = 1                       ! numS=1: spin singlet, numS=2: independent up/down spins 
character(3) :: statename = 'KIV'                        ! label for the solution
integer(dp), parameter:: nenforceC3 = 1                   ! enforce C_3z symmetry
integer(dp), parameter:: nenforceC2 = 0                   ! enforce C_2z symmetry
integer(dp), parameter:: nenforceT = 0                    ! enforce time-reversal (T) symmetry
integer(dp), parameter:: nenforceC2T = 1                  ! enforce C_2zT symmetry
integer(dp), parameter:: nenforceValley = 0               ! enforce U(1)_valley symmetry

! "interaction" parameters
real(dp), parameter :: epsilon = 10.0_dp                   ! dielectric constant
real(dp), parameter :: U = 4.00_dp                       ! Hubbard U in eV
real(dp), parameter :: Delta = 0.000_dp                 ! layer bias in eV
integer(dp), parameter :: numI = 1                      ! number of outer unit cell shells when computing the exchange term
integer(dp), parameter :: numC = 10                      ! number of outer unit cell shells when computing the Hartree term    
integer(dp), parameter :: nscreen = 2                       ! 1/2: single/double metallic gate screening
real(dp), parameter:: xi = 2.0_dp/0.246_dp              ! distance to metallic gates= xi/2, in units of a=2.46 A

! "loop parameters"
real(dp), parameter :: EnergyTolerance=0.0001_dp         ! energies within tolerance are considered degenerate when computing the Fock matrix  
real(dp), parameter :: StepAlternative = 1.0_dp          ! value of step when the ODA algo does not provide it (the energy increases for any value)

! "diagonalization" parameters
integer(dp), parameter :: numb = 5*4                    ! number of bands of the outputs/FockBulk computation/Fermi energy sorting 
integer(dp), parameter :: ncb = numb/2                  ! number of bands above charge neutrality in the outputs (numb-nvb below CN)
integer(dp), parameter :: fphase = 0                    ! 0/1: phase included/not included in the definition of the Bloch states
integer(dp), parameter :: nLower=(ndim - 2*numb + 2*ncb)/2+1    ! Lowest band in outputs/FockBulk computation
integer(dp), parameter :: nUpper=(ndim + 2*ncb)/2               ! Highest band in outputs/FockBulk computation
integer(dp), parameter :: NeutralityPoint=(ndim/2 - nLower)+1   ! Band index  neutrality point

! "output" parameters
character(9) :: dirFock  = 'dataFock/'                                ! folder to output Fock matrix
character(7) :: dir = 'output/'                                       ! folder for all other outputs

! various numerical values
real(dp), parameter :: alpha = 5.853_dp/epsilon*1.0_dp                  ! e^2/(4pi x epsilon0 x epsilon) in units of eV x a    
real(dp), parameter :: alphaH = 5.853_dp/epsilon*1.0_dp                 ! alpha/alphaH is used in the exchange/Hartree term
real(dp), parameter :: pressure=1.0_dp                 ! parameter that mimics hydrostatic pressure. 1.0=no external pressure
real(dp), parameter :: tz = 1.35772_dp/pressure        ! layer separation in units of a
real(dp), parameter :: a0 = 1.0_dp/sqrt(3.0_dp)        ! carbon-carbon distance in units of a
real(dp), parameter :: a1(2) = [ 0.5_dp,sqrt(3.0_dp)*0.5_dp]       ! graphene lattice vectors
real(dp), parameter :: a2(2) = [-0.5_dp,sqrt(3.0_dp)*0.5_dp]       ! graphene lattice vectors
real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)        ! pi=3.141592


! Filename suffixes (built by InitParameters)
character(210) :: parametersOut                          ! parameters of the output
character(210) :: parametersIn                           ! parameters of the input

private :: rstr, istr, cstr, AssignSuffix             ! helpers used to build 'parameters' and 'parametersIn'

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Build 'parameters' and 'parametersIn'.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitParameters()

character(:), allocatable :: s

s = '-'//trim(statename)//'-numS'//istr(numS)//cstr(nenforceC3,nenforceC2,nenforceT,nenforceC2T,nenforceValley)// &
    '-filling'//rstr(nfilling,1,plus=.true.)// &
    '-i'//istr(ntheta)// &
    '-nlayers'//istr(nlayers)// &
    '-TB'//istr(TBFunction)// &
    '-relax'//istr(nrelax)// &
    '-delta'//rstr(Delta*1000.0_dp,1)// &                ! tagged in meV, Delta is in eV
    '-eps'//rstr(epsilon,1)// &
    '-U'//rstr(U,2)// &
    '-screen'//istr(nscreen)// &
    '-xi'//rstr(xi*0.246_dp,1)// &
    '-numI'//istr(numI)// &
    '-numk'//istr(numk)// &
    '-numkPosProc'//istr(numkPostProc)// &
    '-dp'//istr(int(dp,dp))//'.dat'                      ! dp is a default integer, unlike the rest

call AssignSuffix(parametersOut,s,'parameters')

s = '-'//trim(statename)//'-numS'//istr(numS)//cstr(nenforceC3,nenforceC2,nenforceT,nenforceC2T,nenforceValley)// &
    '-filling'//rstr(nfilling,1,plus=.true.)// &
    '-i'//istr(ntheta)// &
    '-nlayers'//istr(nlayers)// &
    '-TB'//istr(TBFunction)// &
    '-relax'//istr(nrelax)// &
    '-delta'//rstr(Delta*1000.0_dp,1)// &                ! tagged in meV, Delta is in eV
    '-eps'//rstr(epsilon,1)// &
    '-U'//rstr(U,2)// &
    '-screen'//istr(nscreen)// &
    '-xi'//rstr(xi*0.246_dp,1)// &
    '-numI'//istr(numI)// &
    '-numk'//istr(numk)// &
    '-dp'//istr(int(dp,dp))//'.dat'                      ! dp is a default integer, unlike the rest

call AssignSuffix(parametersIn,s,'parametersIn')

end subroutine InitParameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! real -> shortest string with nd decimals. A wide Fw.d field keeps the leading zero ('0.0',
!! not '.0' as F0.d would give) and adjustl/trim then removes the padding: no width to overflow.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function rstr(x,nd,plus) result(s)

real(dp), intent(in) :: x
integer, intent(in) :: nd
logical, intent(in), optional :: plus                    ! .true. -> always show the sign, as in '+0.0'
character(:), allocatable :: s
character(40) :: buffer, fmt
logical :: signed

signed = .false.
if(present(plus)) signed = plus

if(signed)then
    write(fmt,'(A,I0,A)') '(SP,F39.',nd,')'
else
    write(fmt,'(A,I0,A)') '(F39.',nd,')'
endif

write(buffer,fmt) x
s = trim(adjustl(buffer))

end function rstr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! integer -> shortest string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function istr(n) result(s)

integer(dp), intent(in) :: n
character(:), allocatable :: s
character(40) :: buffer

write(buffer,'(I0)') n
s = trim(adjustl(buffer))

end function istr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Enforced symmetries as five digits: C3, C2, T, C2T, valley. C2, T and C2T are not independent
!! (C2 x T = C2T), so enforcing any two of them enforces the third and all three are tagged as 1.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function cstr(mC3,mC2,mT,mC2T,mValley) result(s)

integer(dp), intent(in) :: mC3, mC2, mT, mC2T, mValley
character(:), allocatable :: s
integer(dp) :: nC2, nT, nC2T

nC2  = mC2
nT   = mT
nC2T = mC2T

if(nC2+nT+nC2T.EQ.2_dp)then                              ! any two of C2, T, C2T enforce the third
    nC2  = 1_dp
    nT   = 1_dp
    nC2T = 1_dp
endif

s = '-SymConstrain'//istr(mC3)//istr(nC2)//istr(nT)//istr(nC2T)//istr(mValley)

end function cstr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Copies the suffix into its fixed-length variable, stopping instead of silently truncating
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine AssignSuffix(suffix,s,label)

character(*), intent(out) :: suffix
character(*), intent(in) :: s, label

if(len(s).GT.len(suffix))then
    write(*,*) 'ERROR: ',label,' needs ',len(s),' characters, only ',len(suffix),' declared'
    stop 1
endif

suffix = s

end subroutine AssignSuffix


end module Setup
