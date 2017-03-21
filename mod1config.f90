module sysconfig
! module for config of cell cluster
use parameters

contains

! initialize system size and cell position
subroutine initSystem( rCell, rSim, elemMax, pxCell)
    implicit none
    integer, intent(out) :: rCell(:,:), rSim(2), elemMax, pxCell
    integer  :: i, j, k, lCell, lx, ly

    ! calculate cell length in terms of simulation lattices
    lCell = int(sqrt(aCell/pxReal**2))
    if ( lCell == 0 ) then
        lCell = 1
    end if
    ! set up simulation lattice size
    lx = 4 + lfinish       ! x size in terms of cell lengths
    ly = 2                 ! y size in terms of cell lengths
    rSim(1) = lx * lCell
    rSim(2) = ly * lCell

    ! set initial cell location. Start 1 cell radius above x and y boundaries.
    k = 0
    do i = 1, lCell
        do j = 1, lCell
            k = k + 1
            rCell(k,1) = i + (lCell / 2)
            rCell(k,2) = j + (lCell / 2)
        enddo
    enddo
    write(*,*) ' rSim = ', rSim(:)
    write(*,*) ' cell x:', rCell(1,1), rCell(k,1)
    write(*,*) ' cell y:', rCell(1,2), rCell(k,2)
    elemMax = rSim(1) * rSim(2)
    pxCell  = k
end subroutine initSystem


subroutine pickLatticePair( rSim, a, b, rCell, pxCell, elemMax)
    implicit none
    integer, intent(out) :: a(4), b(4)
    integer, intent(in)  :: rSim(2), rCell(:,:), pxCell, elemMax
    integer  :: i
    real(b8) :: r

    ! pick first lattice site and get lattice value
    a(:) = 0
    do i = 1, 2
        call random_number(r)
        a(i) = 1 + floor( real(rSim(i))*r )
    enddo
    do i = 1, pxCell
        if ( a(1) == rCell(i,1) .AND. a(2) == rCell(i,2) ) then
            a(3) = 1
            a(4) = i
            exit
        end if
    enddo

    ! pick second lattice site and get lattice value
    b(:) = 0
    call random_number(r)
    i = 1 + floor(4.0*r)
    call nnGet( i, a, rSim, b)
    if ( b(1) == 0 .OR. b(2) == 0 ) then
        b(3) = 0
    else
        do i = 1, pxCell
            if ( b(1) == rCell(i,1) .AND. b(2) == rCell(i,2) ) then
                b(3) = 1
                b(4) = i
                exit
            end if
        enddo
    end if

end subroutine pickLatticePair


! Output coordinates of the nearest neighbor (nn)
subroutine nnGet( i, x, rSim, nn)
    ! i indicates the nn we are interested in
    ! i = 1 -- nn up
    ! i = 2 -- nn right
    ! i = 3 -- nn down
    ! i = 4 -- nn left
    ! nn = array of nearest neighbor coordinates
    ! rSim = dimensions of simulation space
    ! x = coordinates of point
    implicit none
    integer, intent(in) :: i
    integer, dimension(2), intent(in)  :: rSim, x
    integer, dimension(2), intent(out) :: nn

    if( i == 1 )then
        nn(1) = x(1) + 1
        nn(2) = x(2)
    elseif( i == 2 )then
        nn(1) = x(1)
        nn(2) = x(2) + 1
    elseif( i == 3 )then
        nn(1) = x(1) - 1
        nn(2) = x(2)
    elseif( i == 4 )then
        nn(1) = x(1)
        nn(2) = x(2) - 1
    endif

    if( nn(1) > (rSim(1) + 2) .OR. nn(1) < 1 )then
        nn = [ 0, 0]
    elseif( nn(2) > (rSim(2) + 2) .OR. nn(2) < 1 )then
        nn = [ 0, 0]
    endif

end subroutine nnGet


! delete a pixel from a cell
subroutine delpxCell( rCell, pxCell, px)
    implicit none
    integer, intent(inout) :: rCell(:,:)
    integer, intent(in)    :: pxCell, px(2)
    integer  :: i, j

    do i = 1, pxCell
        if ( px(1) == rCell(i,1) .AND. px(2) == rCell(i,2) ) then
            rCell(i,:) = 0
            do j = i+1, pxCell
                rCell(j-1,:) = rCell(j,:)
            enddo
            exit
        end if
    enddo
end subroutine delpxCell


! Count the number of occupied lattice points
subroutine occupyCount( nl, xcell )
    ! xcell =  array of lattice sites occupied by one cell x(i,:,:)
    ! nl = number of lattice sites contained in xcell
    implicit none
    integer, dimension(:,:), intent(in) :: xcell
    integer, intent(out) :: nl

    ! count how many lattice sites one cell contains
    nl = 1
    do while( xcell(nl,1) /= 0 )
        nl = nl + 1
    enddo
    nl = nl - 1
end subroutine occupyCount


! Checks whether a cell is simply connected or not using flood fill algorithm
recursive subroutine floodFill( node, filled, rSim, xcell)
    ! L = number of lattice sites along one dimension
    ! node = array of lattice site coordinates of node
    ! filled = array of all lattice sites that have been filled by the algorithm
    ! xcell =  array of lattice sites occupied by cell x(i,:,:)
    implicit none
    integer, dimension(2), intent(in) :: node, rSim
    integer, dimension(:,:), intent(inout) :: filled
    integer, dimension(:,:), intent(in) :: xcell
    integer, dimension(2) :: nn
    integer :: i, j, nf, nl

    call occupyCount( nf, filled )
    call occupyCount( nl, xcell)

    do i = 1, nf
        if( node(1) == filled(i,1) .AND. node(2) == filled(i,2) )then
            return
        endif
    enddo

    j = 0
    do i = 1, nl
        if( node(1) == xcell(i,1) .AND. node(2) == xcell(i,2) )then
            j = 1
        endif
    enddo
    if( j /= 1 )then
        return
    endif

    nf = nf + 1
    if( nf > nl )then
        return
    endif

    filled( nf, :) = node

    call nnGet( 1, node, rSim, nn)
    call floodFill( nn, filled, rSim, xcell)

    call nnGet( 2, node, rSim, nn)
    call floodFill( nn, filled, rSim, xcell)

    call nnGet( 3, node, rSim, nn)
    call floodFill( nn, filled, rSim, xcell)

    call nnGet( 4, node, rSim, nn)
    call floodFill( nn, filled, rSim, xcell)

end subroutine floodFill


end module
