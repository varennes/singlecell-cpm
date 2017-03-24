program onecellcpm

! declare modules
use parameters
use sysconfig
use sensing
use probability
use wrtout

! allocate variables
implicit none
integer :: i, j, k, n, nFill
integer :: run, tMC, telem, elemMax

real(b8) :: prob, r, w, ui, uf
real(b8) :: CI, CR, vCell(tMCmax)
real(b8) :: comCell(tMCmax,2), globalSignal
real(b8), allocatable :: localSignal(:)

integer :: a(4), b(4), rSim(2)
integer :: pxCell, pxTmp
integer, allocatable :: rCell(:,:), rTmp(:,:), fill(:,:)

allocate( localSignal( 2*int(aCell/pxReal**2)) )
allocate( rCell( 2*int(aCell/pxReal**2), 2) )
allocate(  rTmp( 2*int(aCell/pxReal**2), 2) )
allocate(  fill( 2*int(aCell/pxReal**2), 2) )

call init_random_seed()

do run = 1, runTotal
    ! run a simulation instance
    comCell = 0.0_b8
    call initSystem( rCell, rSim, elemMax, pxCell)

    ! cell initialization time: cell shape and size relaxes before start of chemotaxis simulation
    do i = 1, elemMax
        call pickLatticePair( rSim, a, b, rCell, pxCell, elemMax)
        if ( a(3) == b(3) .OR. a(1) == 0 .OR. a(2) == 0 .OR. b(1) == 0 .OR. b(2) == 0 ) then
            cycle
        end if
        call getElemStep( a, b, rSim, rCell, pxCell, globalSignal, localSignal)
    enddo

    call getCOM( rCell, comCell(1,:))
    if ( run == 1 ) then
        write(*,*) ' elemMax =', elemMax
        write(*,*) '  pxCell =', pxCell
    end if

    write(106,*) comCell(1,:)
    call wrtCell( rCell, pxCell, 0)

    do tMC = 1, tMCmax
        do telem = 1, elemMax

            ! !! CHECK STATS
            ! if ( tMC == 1 .AND. telem == 1 ) then
            !     write(*,*) 'check 1'
            !     do i = 1, 1000
            !         call getSignal( rCell, pxCell, globalSignal, localSignal)
            !         write(101,*) localSignal(1)
            !         write(102,*) globalSignal
            !     enddo
            ! end if
            ! !! CHECK STATS

            call pickLatticePair( rSim, a, b, rCell, pxCell, elemMax)
            if ( a(3) == b(3) .OR. a(1) == 0 .OR. a(2) == 0 .OR. b(1) == 0 .OR. b(2) == 0 ) then
                ! write(*,*), ' a =', a
                ! write(*,*), ' b =', b
                cycle
            end if
            call getElemStep( a, b, rSim, rCell, pxCell, globalSignal, localSignal)

        enddo ! end of elementary time-step loop
        call getCOM( rCell, comCell(tMC,:))

        write(106,*) comCell(tMC,:)
        call wrtCell( rCell, pxCell, tMC)

    enddo ! end of Monte Carlo time-step loop

    call getCellSpeed( tMC-1, 3, comCell, vCell)
    call getChemotaxMetric( tMC-1, comCell, CI, CR)
    write(121,*) CI, CR, run

    i = 1
    do while ( vCell(i) /= 0.0_b8 )
        write(122,*) vCell(i), run
        i = i + 1
    enddo

enddo

deallocate( localSignal )
deallocate( rCell )
deallocate( rTmp  )

end program

! subroutine: initialize RANDOM_SEED
subroutine init_random_seed()
    integer :: values(1:8), k
    integer, dimension(:), allocatable :: seed
    real(8) :: r

    call date_and_time(values=values)

    call random_seed(size=k)
    allocate(seed(1:k))
    seed(:) = values(8)
    call random_seed(put=seed)
end subroutine init_random_seed
