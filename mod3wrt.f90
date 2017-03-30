module wrtout
! module for output of simulation information
use parameters

contains

! outputs x to fort.105
subroutine wrtCell( rCell, pxCell, t)
    implicit none
    integer, intent(in) :: rCell(:,:), pxCell, t
    integer :: i, j

    do i = 1, pxCell
        write(105,*) rCell(i,1), rCell(i,2), t
    enddo
end subroutine wrtCell


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


subroutine wrtInstSpeed( vCell, dt, run)
    implicit none
    real(b8), intent(in) :: vCell(:)
    integer,  intent(in) :: dt, run
    integer :: i, t
    character(len=1024) :: filename

    write (filename,"(A10)") 'v_inst.dat'
    open( 14, file=filename)

    i = 1
    t = i + dt
    do while ( vCell(i) /= 0.0_b8 )
        write(14,"(E16.8)", advance="no") vCell(i)
        write(14,"(I7)", advance="no") t
        write(14,"(I7)", advance="no") run
        write(14,*) ''
        i = i + 1
        t = t + dt
    enddo
end subroutine wrtInstSpeed


end module
