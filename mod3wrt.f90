module wrtout
! module for output of simulation information
use parameters
use probability

contains

! outputs cell pixel locations to fort.105
subroutine wrtCell( rCell, pxCell, t)
    implicit none
    integer,  intent(in) :: rCell(:,:), pxCell, t
    integer :: i
    do i = 1, pxCell
        write(105,*) rCell(i,1), rCell(i,2), t
    enddo
end subroutine wrtCell


! outputs center-of-mass (COM) to fort.106
subroutine wrtCOM( comCell, t, run)
    implicit none
    real(b8), intent(in) :: comCell(2)
    integer,  intent(in) :: t, run
    character(len=1024)  :: filename
    write (filename,"(A3,I0.3,A4)") 'com', run, '.dat'
    open( 22, file=filename)

    write(22,"(E16.8)", advance="no") comCell(1)
    write(22,"(E16.8)", advance="no") comCell(2)
    write(22,"(I7)", advance="no") t
    write(22,*) ''
end subroutine wrtCOM


! output x, y displacement and squared displacement
subroutine wrtDisplacement( comCell, com0, t, run)
    implicit none
    real(b8), intent(in) :: comCell(2), com0(2)
    integer,  intent(in) :: run, t
    character(len=1024)  :: filename
    real(b8) :: r2
    r2 = (comCell(1) - com0(1))**2 + (comCell(2) - com0(2))**2
    write (filename,"(A11)") 'dsplcmt.dat'
    open( 21, file=filename)

    write(21,"(E16.8)", advance="no") comCell(1) - com0(1)
    write(21,"(E16.8)", advance="no") comCell(2) - com0(2)
    write(21,"(E16.8)", advance="no") r2
    write(21,"(I7)",    advance="no") t
    write(21,"(I7)",    advance="no") run
    write(21,*) ''
end subroutine wrtDisplacement


subroutine wrtChemotaxMetric( CI, CR, run)
    implicit none
    real(b8), intent(in) :: CI, CR
    integer,  intent(in) :: run
    character(len=1024) :: filename

    write (filename, "(A6)") 'ci.dat'
    open( 12, file=filename)
    write (filename, "(A6)") 'cr.dat'
    open( 13, file=filename)

    write(12,"(E16.8)", advance="no") CI
    write(12,"(I7)", advance="no") run
    write(12,*) ''

    write(13,"(E16.8)", advance="no") CR
    write(13,"(I7)", advance="no") run
    write(13,*) ''
end subroutine wrtChemotaxMetric


subroutine wrtInstSpeed( vCell, dt, run, iv)
    implicit none
    real(b8), intent(in) :: vCell(:)
    integer,  intent(in) :: dt, run, iv
    integer :: i, t
    character(len=1024) :: filename

    write (filename,"(A10)") 'v_inst.dat'
    open( 14, file=filename)

    t = 1 + dt
    do i = 1, iv
        write(14,"(E16.8)", advance="no") vCell(i)
        write(14,"(I7)", advance="no") t
        write(14,"(I7)", advance="no") run
        write(14,*) ''
        t = t + dt
    enddo
end subroutine wrtInstSpeed


subroutine wrtMeanSpeed( vCell, run, iv)
    implicit none
    real(b8), intent(in) :: vCell(:)
    integer,  intent(in) :: run, iv
    real(b8) :: mean
    character(len=1024) :: filename

    write (filename,"(A10)") 'v_mean.dat'
    open( 15, file=filename)

    mean = sum(vCell(:))
    mean = mean / float(iv)

    write(15,"(E16.8)", advance="no") mean
    write(15,"(I7)", advance="no") run
    write(15,*) ''
end subroutine wrtMeanSpeed


subroutine wrtPlrVec( plrVec, t, run)
    implicit none
    real(b8), intent(in) :: plrVec(2)
    integer,  intent(in) :: t, run
    character(len=1024)  :: filename
    write (filename,"(A3,I0.3,A4)") 'plr', run, '.dat'
    open( 33, file=filename)

    write(33,"(E16.8)", advance="no") plrVec(1)
    write(33,"(E16.8)", advance="no") plrVec(2)
    write(33,"(I7)", advance="no") t
    write(33,*) ''
end subroutine wrtPlrVec


! write out cell aread and perimeter in simulation units
subroutine wrtAspectRatio( pxCell, rCell, run)
    implicit none
    integer, intent(in)  :: rCell(:,:), pxCell, run
    integer :: perim
    character(len=1024) :: filename

    call perimCheck( rCell, pxCell, perim)

    write (filename,"(A10)") 'aspect.dat'
    open( 16, file=filename)

    write(16,"(I7)", advance="no") pxCell
    write(16,"(I7)", advance="no") perim
    write(16,"(I7)", advance="no") run
    write(16,*) ''
end subroutine wrtAspectRatio


end module
