module sysconfig
! module for config of cell cluster
use parameters

contains

! initialize system size and cell position
subroutine initSystem( rCell, rSim, elemMax, pxCell, pCell)
    implicit none
    real(b8), intent(out) :: pCell
    integer,  intent(out) :: rCell(:,:), rSim(2), elemMax, pxCell
    integer  :: i, j, k, lCell, lx, ly

    ! calculate cell length in terms of simulation lattices
    lCell = int(dsqrt(aCell/pxReal**2))
    if ( lCell == 0 ) then
        lCell = 1
    end if
    ! set up simulation lattice size
    lx = 4 + lfinish       ! x size in terms of cell lengths
    ly = 7                 ! y size in terms of cell lengths
    rSim(1) = lx * lCell
    rSim(2) = ly * lCell

    ! set initial cell location. Start 1 cell radius above x and y boundaries.
    k = 0
    do i = 1, lCell
        do j = 1, lCell
            k = k + 1
            rCell(k,1) = i + 2*lCell
            rCell(k,2) = j + ( rSim(2) - lCell) / 2
        enddo
    enddo
    pCell = 2.0_b8 * dsqrt( pi * aCell)
    elemMax = rSim(1) * rSim(2)
    pxCell  = k

end subroutine initSystem


! calculate cell centroid (aka COM)
subroutine getCOM( rCell, com)
    implicit none
    integer,  intent(in)  :: rCell(:,:)
    real(b8), intent(out) :: com(2)
    integer  :: i
    real(b8) :: x, y

    x = 0.0_b8
    y = 0.0_b8
    i = 1
    do while( rCell(i,1) /= 0 )
        x = x + real(rCell(i,1))
        y = y + real(rCell(i,2))
        i = i + 1
        if ( i > 4*int(aCell/pxReal**2) ) then
            exit
        end if
    enddo
    i = i - 1
    com(1) = x / real(i)
    com(2) = y / real(i)
end subroutine getCOM


subroutine getChemotaxMetric( tf, dt, comCell, CI, CR)
    implicit none
    integer,  intent(in)  :: tf, dt
    real(b8), intent(in)  :: comCell(:,:)
    real(b8), intent(out) :: CI, CR
    real(b8) :: displacement, distance
    integer  :: t

    displacement = dsqrt( (comCell(tf,1) - comCell(1,1))**2 + (comCell(tf,2) - comCell(1,2))**2 )
    distance = 0.0_b8
    do t = 1+dt, tf, dt
        distance = distance + dsqrt( (comCell(t,1) - comCell(t-dt,1))**2 + (comCell(t,2) - comCell(t-dt,2))**2 )
    enddo

    if ( distance > 1E-10 ) then
        CR = displacement / distance
    else
        CR = 0.0_b8
    end if

    if ( displacement > 1E-10 ) then
        CI = (comCell(tf,1) - comCell(1,1)) / displacement
    else
        CI = 0.0_b8
    end if
end subroutine getChemotaxMetric


subroutine getCellSpeed( tf, dt, comCell, vCell, iv)
    implicit none
    integer,  intent(in)  :: tf, dt
    integer,  intent(out) :: iv
    real(b8), intent(in)  :: comCell(:,:)
    real(b8), intent(out) :: vCell(:)
    real(b8) :: d
    integer  :: i, t

    vCell = 0.0_b8
    i = 0
    do t = 1+dt, tf, dt
        i = i + 1
        d = dsqrt( (comCell(t,1) - comCell(t-dt,1))**2 + (comCell(t,2) - comCell(t-dt,2))**2 )
        vCell(i) = d / real(dt)
    enddo
    iv = i
    ! write(*,*) t, tf, iv, vCell(iv), vCell(iv-1), vCell(iv+1)
end subroutine getCellSpeed


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


! delete a pixel from a cell
subroutine delpxCell( rCell, pxCell, px)
    implicit none
    integer, intent(inout) :: rCell(:,:)
    integer, intent(in)    :: pxCell, px(2)
    integer  :: i, j

    do i = 1, pxCell
        if ( px(1) == rCell(i,1) .AND. px(2) == rCell(i,2) ) then
            ! write(*,*) '   deleted, i=', i, '|', rCell(i,1:2)
            rCell(i,:) = 0
            do j = i+1, pxCell
                rCell(j-1,:) = rCell(j,:)
            enddo
            rCell(pxCell,:) = 0
            exit
        end if
    enddo
end subroutine delpxCell


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
