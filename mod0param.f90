module parameters

! b8 will be used to define reals with 14 digits
integer,  parameter :: b8 = selected_real_kind(14)

!!! REAL UNIT PARAMETERS !!!
real(b8), parameter :: alpha =   1.0_b8   ! in energy
real(b8), parameter :: lArea =   1.0_b8   ! in energy
real(b8), parameter :: w0    =   1.0_b8  ! in energy
real(b8), parameter :: aCell = 100.0_b8   ! in microns^2
real(b8), parameter :: c0    =  10.0_b8   ! in nanoMolars
real(b8), parameter :: g     =   1.0_b8   ! in nanoMolars / micron

!!! SIMULATION PARAMETERS
real(b8), parameter :: pxReal   =  1.0_b8 ! length of one lattice site in microns
integer,  parameter :: lfinish  =  5      ! finish line in terms of cell lengths
integer,  parameter :: runTotal =  1      ! total number of runs
integer,  parameter :: tMCmax   = 1      ! max number of MC time-steps

end module

!!! PARAMETER DESCRIPTIONS !!!
!
! alpha = cell-ECM contact energy
! lArea = Area fluctuation energy cost
!    w0 = Work proportionality factor
! aCell = Cell resting area
!    c0 = mean chemical concentration at x=0
!     g = mean gradient
!
!!! PARAMETER DESCRIPTIONS !!!
