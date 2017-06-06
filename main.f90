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
real(b8) :: pCell, plrVec(2), globalSignal
real(b8) :: comCell(tMCmax,2), deltaCOM(2), comOld(2), comNew(2)
real(b8), allocatable :: localSignal(:)

integer :: a(4), b(4), rSim(2)
integer :: pxCell, pxTmp
integer, allocatable :: rCell(:,:), rTmp(:,:), fill(:,:)

character(len=1024) :: filename

allocate( localSignal( 4*int(aCell/pxReal**2)) )
allocate( rCell( 4*int(aCell/pxReal**2), 2) )
allocate(  rTmp( 4*int(aCell/pxReal**2), 2) )
allocate(  fill( 4*int(aCell/pxReal**2), 2) )

call init_random_seed()

do run = 1, runTotal
    ! run a simulation instance
    comCell = 0.0_b8
    call initSystem( rCell, rSim, elemMax, pxCell, pCell)

    if ( run == 1 ) then
        write(*,*) ' rSim = ', rSim(:)
        write(*,*) ' cell x:', rCell(1,1), rCell(pxCell,1)
        write(*,*) ' cell y:', rCell(1,2), rCell(pxCell,2)
        write(*,*) ' pxCell:   ',   pxCell, '| elemMax =', elemMax
        write(*,*) ' runTotal =', runTotal, '| tMCmax =', tMCmax
        write(*,*)
    end if

    ! cell initialization time: cell shape and size relaxes before start of chemotaxis simulation
    plrVec = 0.0_b8
    comOld = 0.0_b8
    comNew = 0.0_b8
    do i = 1, 4*elemMax
        call pickLatticePair( rSim, a, b, rCell, pxCell, elemMax)
        if ( a(3) == b(3) .OR. a(1) == 0 .OR. a(2) == 0 .OR. b(1) == 0 .OR. b(2) == 0 ) then
            cycle
        end if
        call getItlStep( a, b, rSim, rCell, pxCell, pCell)
        call getCOM( rCell, comNew)
        deltaCOM = comNew - comOld
        call getVectorUpdate( plrVec, rCell, pxCell, comNew, deltaCOM, globalSignal, localSignal, rSim)
        comOld = comNew
    enddo

    call getCOM( rCell, comCell(1,:))

    ! write(107,*) plrVec, 1
    ! call wrtCell( rCell, comCell(1,:), pxCell, 1)

    do tMC = 2, tMCmax
        call getCOM( rCell, comCell(tMC,:))
        deltaCOM = comCell(tMC,:) - comCell(tMC-1,:)
        call getVectorUpdate( plrVec, rCell, pxCell, comCell(tMC,:), deltaCOM, globalSignal, localSignal, rSim)
        do telem = 1, elemMax

            call pickLatticePair( rSim, a, b, rCell, pxCell, elemMax)
            if ( a(3) == b(3) .OR. a(1) == 0 .OR. a(2) == 0 .OR. b(1) == 0 .OR. b(2) == 0 ) then
                cycle
            end if
            call getVectorStep( a, b, rSim, rCell, pxCell, pCell, plrVec)

        enddo ! end of elementary time-step loop

        if ( mod(tMC,3) == 0 ) then
            ! write(107,*) plrVec, tMC
            call wrtCell( rCell, comCell(tMC,:), pxCell, tMC, run)
        end if

        ! check if finish line hit
        if ( comCell(tMC,1) > (real(rSim(1))-sqrt(aCell/pxReal**2)) ) then
            write(*,*) '  finish line hit'
            exit
        end if
    enddo ! end of Monte Carlo time-step loop

    dt = 3
    call getCellSpeed( tMC-1, dt, comCell, vCell, iv)
    call getChemotaxMetric( tMC-1, comCell, CI, CR)

    ! call wrtDisplacement( comCell(tMC-1,:), comCell(1,:), tMC-1, run)
    call wrtChemotaxMetric( CI, CR, run)
    call wrtMeanSpeed( vCell, run, iv)
    ! call wrtInstSpeed( vCell, dt, run, iv)

    write(*,"(A14)", advance="no") 'complete run #'
    write(*,"(I5)") run

enddo

close(12)
close(13)
close(14)
close(15)
close(16)
close(21)
close(22)

! deallocate( localSignal )
! deallocate( rCell )
! deallocate( rTmp  )

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
