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

integer :: a(4), b(4), rSim(2,2)
integer :: pxCell, pxTmp
integer, allocatable :: rCell(:,:), rTmp(:,:), fill(:,:)

character(len=1024) :: filename

real(b8) :: cpuT0, cpuT1

call cpu_time(cpuT0)

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
        write(*,*) ' rSim(1) = ', rSim(1,:)
        write(*,*) ' rSim(2) = ', rSim(2,:)
        write(*,*) ' cell x:', rCell(1,1), rCell(pxCell,1)
        write(*,*) ' cell y:', rCell(1,2), rCell(pxCell,2)
        write(*,*) ' pxCell:   ',   pxCell, '| elemMax =', elemMax
        write(*,*) ' runTotal =', runTotal, '| tMCmax =', tMCmax
        write(*,*)
    end if

    ! initialize polarization and COM values
    plrVec = 0.0_b8
    comOld = 0.0_b8
    comNew = 0.0_b8
    ! initialize cell; cell shape relaxes
    do i = 1, 10
        ! update elemMax
        do k = 1, 2
            rSim(k,1) = minval( rCell(1:pxCell,k)) - 1
            rSim(k,2) = maxval( rCell(1:pxCell,k)) + 1
        enddo
        elemMax = (rSim(1,2) - rSim(1,1) + 1) * (rSim(2,2) - rSim(2,1) + 1)

        do j = 1, elemMax
            call pickLatticePair( rSim, a, b, rCell, pxCell, elemMax)
            if ( a(3) == b(3) .OR. a(1) == 0 .OR. a(2) == 0 .OR. b(1) == 0 .OR. b(2) == 0 ) then
                cycle
            end if
            call getItlStep( a, b, rCell, pxCell, pCell)
            call getCOM( rCell, comNew)
            deltaCOM = comNew - comOld
        enddo
        ! call getVectorUpdate( plrVec, rCell, pxCell, comNew, deltaCOM, globalSignal, localSignal)
        comOld = comNew
    enddo
    plrVec = 0.0_b8

    ! get initial cell COM
    call getCOM( rCell, comCell(1,:))

    ! ! output initial cell information
    ! call wrtCell( rCell, pxCell, 1)
    ! call wrtCOM( comCell(1,:), 1, run)
    ! call wrtPlrVec( plrVec, 1, run)

    do tMC = 2, tMCmax
        ! update location of COM and change in COM
        call getCOM( rCell, comCell(tMC,:))
        deltaCOM = comCell(tMC,:) - comCell(tMC-1,:)

        ! update polarization vector
        call getVectorUpdate( plrVec, rCell, pxCell, comCell(tMC,:), deltaCOM, globalSignal, localSignal)

        ! update elemMax
        do i = 1, 2
            rSim(i,1) = minval( rCell(1:pxCell,i)) - 1
            rSim(i,2) = maxval( rCell(1:pxCell,i)) + 1
        enddo
        elemMax = (rSim(1,2) - rSim(1,1) + 1) * (rSim(2,2) - rSim(2,1) + 1)

        ! elementary time-step loop
        do telem = 1, elemMax
            call pickLatticePair( rSim, a, b, rCell, pxCell, elemMax)
            if ( a(3) == b(3) ) then
                cycle
            end if
            call getVectorStep( a, b, rCell, pxCell, pCell, plrVec)
        enddo

        ! ! output cell information
        ! if ( mod(tMC, 10) == 0 ) then
        !     call wrtCell( rCell, pxCell, tMC)
        !     call wrtCOM( comCell(tMC,:), tMC, run)
        !     call wrtPlrVec( plrVec, tMC, run)
        ! end if

    enddo ! end of Monte Carlo time-step loop

    dt = 10
    call getCellSpeed( tMC-1, dt, comCell, vCell, iv)
    call getChemotaxMetric( tMC-1, dt, comCell, CI, CR)

    call wrtChemotaxMetric( CI, CR, run)
    call wrtMeanSpeed( vCell, run, iv)
    ! call wrtDisplacement( comCell(tMC-1,:), comCell(1,:), tMC-1, run)
    ! call wrtInstSpeed( vCell, dt, run, iv)

    write(*,"(A14)", advance="no") 'complete run #'
    write(*,"(I5)") run

enddo

call cpu_time(cpuT1)

write(*,*) '  SIMULATION COMPLETE'
write(*,*) ''
write(*,*) '  run time (min) =', ( cpuT1 - cpuT0) / 60.0_b8

close(12)
close(13)
close(14)
close(15)
close(16)
close(21)
close(22)


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
