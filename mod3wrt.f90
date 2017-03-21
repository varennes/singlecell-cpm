module wrtout
! module for output of simulation information

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

end module
