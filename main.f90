program onecellcpm

! declare modules
use parameters
use sysconfig
use sensing
use probability
use wrtout

! allocate variables
implicit none
integer :: i, iv, j, k, n, nFill
integer :: run, dt, tMC, telem, elemMax

real(b8) :: prob, r, w, ui, uf
real(b8) :: CI, CR, vCell(tMCmax)
real(b8) :: comCell(tMCmax,2), pCell, globalSignal
real(b8), allocatable :: localSignal(:)

integer :: a(4), b(4), rSim(2)
integer :: pxCell, pxTmp
integer, allocatable :: rCell(:,:), rTmp(:,:), fill(:,:)

character(len=1024) :: filename

allocate( localSignal( 2*int(aCell/pxReal**2)) )
allocate( rCell( 2*int(aCell/pxReal**2), 2) )
allocate(  rTmp( 2*int(aCell/pxReal**2), 2) )
allocate(  fill( 2*int(aCell/pxReal**2), 2) )

call init_random_seed()

do run = 1, runTotal
    ! run a simulation instance
    comCell = 0.0_b8
    call initSystem( rCell, rSim, elemMax, pxCell, pCell)

    if ( run == 1 ) then
        write(*,*) ' rSim = ', rSim(:)
        write(*,*) ' cell x:', rCell(1,1), rCell(pxCell,1)
        write(*,*) ' cell y:', rCell(1,2), rCell(pxCell,2)
        write(*,*) ' pxCell:', pxCell
        write(*,*) ' elemMax =', elemMax
        write(*,*)
    end if

    ! cell initialization time: cell shape and size relaxes before start of chemotaxis simulation
    do i = 1, 4*elemMax
        call pickLatticePair( rSim, a, b, rCell, pxCell, elemMax)
        if ( a(3) == b(3) .OR. a(1) == 0 .OR. a(2) == 0 .OR. b(1) == 0 .OR. b(2) == 0 ) then
            cycle
        end if
        call getItlStep( a, b, rSim, rCell, pxCell, pCell)

    enddo

    call getCOM( rCell, comCell(1,:))

    call wrtCell( rCell, comCell(1,:), pxCell, 0)

    do tMC = 1, tMCmax
        do telem = 1, elemMax

            call pickLatticePair( rSim, a, b, rCell, pxCell, elemMax)
            if ( a(3) == b(3) .OR. a(1) == 0 .OR. a(2) == 0 .OR. b(1) == 0 .OR. b(2) == 0 ) then
                cycle
            end if
            call getElemStep( a, b, rSim, rCell, pxCell, pCell, globalSignal, localSignal)

        enddo ! end of elementary time-step loop
        call getCOM( rCell, comCell(tMC,:))

        if ( mod(tMC,10) == 0 ) then
            call wrtCell( rCell, comCell(tMC,:), pxCell, tMC)
        end if

        ! check if finish line hit
        if ( comCell(tMC,1) > (real(rSim(1))-sqrt(aCell/pxReal**2)) ) then
            write(*,*) '  finish line hit'
            exit
        end if
    enddo ! end of Monte Carlo time-step loop

    dt = 5
    call getCellSpeed( tMC-1, dt, comCell, vCell, iv)
    call getChemotaxMetric( tMC-1, comCell, CI, CR)

    ! call wrtChemotaxMetric( CI, CR, run)
    ! call wrtMeanSpeed( vCell, run, iv)
    ! call wrtInstSpeed( vCell, dt, run, iv)

    write(*,"(A14)", advance="no") 'complete run #'
    write(*,"(I5)") run

enddo

close(12)
close(13)
close(14)
close(15)

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
