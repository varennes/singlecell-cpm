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
real(b8) :: globalSignal
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
    call initSystem( rCell, rSim, elemMax, pxCell)
    write(*,*) ' elemMax =', elemMax
    write(*,*) '  pxCell =', pxCell

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
            if ( a(3) == b(3) ) then
                ! write(*,*), ' a =', a
                ! write(*,*), ' b =', b
                cycle
            end if
            call getSignal( rCell, pxCell, globalSignal, localSignal)

            if ( a(3) == 1 ) then
                ! cell is adding a pixel
                rTmp = 0
                rTmp = rCell
                pxTmp = pxCell + 1
                rTmp(pxTmp,1:2) = b(1:2)

                w = getWork( globalSignal, localSignal(a(4)), a(3))

                ui = getEnergy( rSim, rCell, pxCell)
                uf = getEnergy( rSim, rTmp,  pxTmp)

                prob = getProb( uf, ui, w)
            else
                ! cell is removing a pixel
                fill = 0
                rTmp = 0
                rTmp = rCell
                pxTmp = pxCell - 1
                call delpxCell( rTmp, pxTmp+1, b(1:2))
                ! check if cell pixels are simply connected
                call floodFill( rTmp(1,1:2), fill, rSim, rTmp(:,:))
                call occupyCount( nFill, fill )
                if ( nFill /= pxTmp .OR. pxTmp == 0 ) then
                    ! cell is not simply connected
                    ! write(*,*) 'nF=', nFill, 'pxTmp=', pxTmp
                    prob = 0.0_b8
                else
                    w = getWork( globalSignal, localSignal(b(4)), a(3))

                    ui = getEnergy( rSim, rCell, pxCell)
                    uf = getEnergy( rSim, rTmp,  pxTmp)

                    prob = getProb( uf, ui, w)
                end if

            end if

            ! write(*,*) ' prob, w = ', prob, w

            call random_number(r)
            if ( r < prob ) then
                 rCell =  rTmp
                pxCell = pxTmp
                ! write(*,*) 'success: a, b', a, b
            end if


        enddo ! end of elementary time-step loop
    call wrtCell( rCell, pxCell, tMC)

    enddo ! end of Monte Carlo time-step loop

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
