
module Setup

implicit none

integer, parameter :: dp = 8                                                               ! dp=4/8: single/double precision
integer(dp), parameter :: nlayers = 2                                                      ! number of layers
integer(dp), parameter :: ntheta = 9                                                       ! twist angle = acos(1 - 1/(6ntheta^2 + 6ntheta + 2))
integer(dp), parameter :: RotateLayers(nlayers) = [-1,+1]                                  ! Positive/negative rotation of layers 
integer(dp), parameter :: ndim = nlayers*2*(ntheta**2+(ntheta+1)**2+(ntheta+1)*ntheta)     ! number of atoms in unt cell
integer(dp), parameter :: numk = 6                       ! number of grid points in the BZ  Nk = numk * numk
integer(dp), parameter :: nrelax = 0                     ! 0/1: unrelaxed/relaxed lattice  
integer(dp), parameter :: numS = 1                       ! numS=1: spin singlet, numS=2: independent up/down spins 
real(dp), parameter :: nfilling = 0.0_dp                 ! filling in units of e/unit cell. -4/+4:empty/full flat bands, 0:charge neutrality
character(3) :: statename = 'SSS'                        ! label for the solution
real(dp), parameter :: U = 4.00_dp                       ! Hubbard U in eV
real(dp), parameter :: epsilon = 10.0                   ! dielectric constant
real(dp), parameter :: Delta = 0.000_dp                 ! layer bias in eV
integer(dp), parameter :: numI = 1                      ! number of outer unit cell shells when computing the exchange term
integer(dp), parameter :: numC = 3                      ! number of outer unit cell shells when computing the Hartree term    
integer(dp), parameter :: nscreen = 2                       ! 1/2: single/double metallic gate screening
real(dp), parameter:: xi = 10.0_dp/0.246_dp              ! distance to metallic gates= xi/2, in units of a=2.46 A
integer(dp), parameter :: numb = minval([20_dp,ndim])   ! number of bands of the outputs/FockBulk computation/Fermi energy sorting 
integer(dp), parameter :: ncb = numb/2                  ! number of bands above charge neutrality in the outputs (numb-nvb below CN)
integer(dp), parameter :: fphase = 0                    ! 0/1: phase included/not included in the definition of the Bloch states

integer(dp), parameter:: nenforceC3 = 0                   ! enforce C_3z symmetry
integer(dp), parameter:: nenforceC2 = 0                   ! enforce C_2z symmetry
integer(dp), parameter:: nenforceT = 0                    ! enforce time-reversal (T) symmetry
integer(dp), parameter:: nenforceC2T = 0                  ! enforce C_2zT symmetry
integer(dp), parameter:: nenforceValley = 0               ! enforce U(1)_valley symmetry

integer(dp), parameter :: nPrintOutputs = 1              ! 0/1: do not print/print outputs
integer(dp), parameter :: nWriteFock = 0                 ! 0/1: do not write/write Fock matrix as an output

integer(dp), parameter :: nRead = 0                      ! 0/1: compute/read inital guess for the state 
integer(dp), parameter :: dpIn = 8                       ! precision of initial guess if nRead=1
integer(dp), parameter :: numkIn = 6                       ! numk of initial guess
integer(dp), parameter :: nrelaxIn = 0                     ! nrelax of initial guess
integer(dp), parameter :: numSIn = 1                       ! numS of initial guess 
real(dp), parameter :: nfillingIn = 0.0_dp                 ! filling of initial guess
character(3) :: statenameIn = 'SSS'                        ! label for the solution of initial guess
real(dp), parameter :: UIn = 4.00_dp                       ! Hubbard U in eV of initial guess
real(dp), parameter :: epsilonIn = 10.0                    ! dielectric constant of initial guess
real(dp), parameter :: DeltaIn = 0.000_dp                 ! layer bias in eV of initial guess
integer(dp), parameter :: numIIn = 1                      ! numI of initial guess
integer(dp), parameter :: numCIn = 3                      ! numc of initial guess
integer(dp), parameter :: nscreenIn = 2                       ! scr of initial guess
real(dp), parameter:: xiIn = 2.0_dp/0.246_dp              ! xi of initial guess

real(dp), parameter :: EnergyTolerance=0.0001_dp         ! energies within tolerance are considered degenerate when computing the Fock matrix  
real(dp), parameter :: StepAlternative = 1.0_dp          ! value of step when the ODA algo does not provide it (the energy increases for any value)
integer(dp), parameter :: itmax = 2                      ! maximal number of iterations of the self-consistency loop  

character(9) :: dirFock  = 'dataFock/'                                ! folder to output Fock matrix
character(7) :: dir = 'output/'                                       ! folder for all other outputs

real(dp), parameter :: alpha = 5.853_dp/epsilon*1.0_dp                  ! e^2/(4pi x epsilon0 x epsilon) in units of eV x a    
real(dp), parameter :: alphaH = 5.853_dp/epsilon*1.0_dp                 ! alpha/alphaH is used in the exchange/Hartree term

real(dp), parameter :: pressure=1.0_dp                 ! parameter that mimics hydrostatic pressure. 1.0=no external pressure
real(dp), parameter :: a0 = 1.0_dp/sqrt(3.0_dp)        ! carbon-carbon distance in units of a
real(dp), parameter :: r0 = 0.184_dp                   ! tight-binding decay constant
real(dp), parameter :: d0 = 1.35772_dp/pressure        ! tight-binding decay constant
real(dp), parameter :: tz = 1.35772_dp/pressure        ! layer separation in units of a

integer(dp), parameter :: nLower=(ndim - 2*numb + 2*ncb)/2+1    ! Lowest band in outputs/FockBulk computation
integer(dp), parameter :: nUpper=(ndim + 2*ncb)/2               ! Highest band in outputs/FockBulk computation
integer(dp), parameter :: NeutralityPoint=(ndim/2 - nLower)+1   ! Band index  neutrality point

real(dp), parameter :: a1(2) = [ 0.5_dp,sqrt(3.0_dp)*0.5_dp]       ! graphene lattice vectors
real(dp), parameter :: a2(2) = [-0.5_dp,sqrt(3.0_dp)*0.5_dp]       ! graphene lattice vectors

real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)        ! pi=3.141592

end module Setup
