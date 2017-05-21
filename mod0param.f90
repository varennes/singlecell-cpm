module parameters

! b8 will be used to define reals with 14 digits
integer,  parameter :: b8 = selected_real_kind(14)

!!! REAL UNIT PARAMETERS !!!
real(b8), parameter :: alpha =  0.5_b8
real(b8), parameter :: lArea =  0.3_b8
real(b8), parameter :: lPerim =  0.01_b8
real(b8), parameter :: w0 =  0.5_b8
real(b8), parameter :: aCell =  400.0_b8
real(b8), parameter :: c0 =  0.0_b8
real(b8), parameter :: g =  1.0_b8

!!! SIMULATION PARAMETERS
real(b8), parameter :: pxReal =  3.0_b8
integer,  parameter :: lfinish =  3
integer,  parameter :: runTotal =  1
integer,  parameter :: tMCmax =  500

real(b8), parameter :: pi = 3.1415927_b8  ! pi

end module

!!! PARAMETER DESCRIPTIONS !!!
!
! alpha  = cell-ECM contact energy in units of energy
! lArea  = Area fluctuation energy cost in units of energy
! lPerim = Perimeter fluctuation energy cost in units of energy
!    w0  = Work proportionality factor in units of energy
! aCell  = Cell resting area in units of microns^2
!    c0  = mean chemical concentration at x=0 in units of nanoMolars
!     g  = mean gradient in units of nanoMolars / micron
!
!   pxReal = length of one lattice site in microns
!  lfinish = finish line in terms of cell lengths
! runTotal = total number of runs
!   tMCmax = max number of MC time-steps
!
!!! PARAMETER DESCRIPTIONS !!!
