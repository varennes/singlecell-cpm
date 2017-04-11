module probability
! module for calculating probability
use parameters
use sysconfig
use sensing

contains

real(b8) function getWork( global, local, sigma)
    implicit none
    real(b8), intent(in) :: global, local
    integer,  intent(in) :: sigma

    getWork = w0 * (local-global) * (-1.0_b8)**(1-sigma)
end function getWork


! evaluate energy of the configuration
real(b8) function getEnergy( rSim, rCell, pxCell, pCell)
    implicit none
    real(b8), intent(in) :: pCell
    integer,  intent(in) :: rSim(2), rCell(:,:), pxCell
    integer :: nn(2)
    integer :: i, j, perim
    real(b8) :: contact, areaCost, perimCost

    ! get cell perimeter for contact energy and perimeter energy
    call perimCheck( rCell, pxCell, rSim, perim)
    contact = alpha * pxReal**2 * real(perim)
    perimCost = lPerim * pxReal * ( real(perim) - (pCell/pxReal))**2

    ! energy contribution due to area
    areaCost = lArea * pxReal**4 * ( real(pxCell) - (aCell/pxReal**2))**2
    getEnergy = contact + areaCost + perimCost

end function getEnergy


! calculate the probability of step being executed
real(b8) function getProb( uNew, uOld, w)
    implicit none
    real(b8), intent(in) :: uNew, uOld, w
    real(b8) :: du

    ! calculate the change in the goal function
    du = uNew - uOld
    ! calculate probability
    getProb = exp( min(0.0, w - du ) )

end function getProb


! energy contribution due to J
real(b8) function jCheck( lCell, lnn, rCell, pxCell)
    ! lCell = lattice site from cell list
    ! lnn   = lattice site nn of lCell
    implicit none
    integer, intent(in) :: lCell(2), lnn(2), rCell(:,:), pxCell
    integer :: i, nnSigma

    nnSigma = 0
    do i = 1, pxCell
        if ( lnn(1) == rCell(i,1) .AND. lnn(2) == rCell(i,2) ) then
            nnSigma = 1
            exit
        end if
    enddo

    if ( nnSigma == 0 ) then
        jCheck = alpha
    else
        jCheck = 0.0_b8
    endif
end function jCheck


! subroutine of all the steps neccasary for one elementary time-step
subroutine getElemStep( a, b, rSim, rCell, pxCell, pCell, globalSignal, localSignal)
    implicit none
    integer, intent(in)     :: a(4), b(4), rSim(2)
    integer, intent(inout)  :: rCell(:,:), pxCell
    real(b8), intent(inout) :: globalSignal, localSignal(:), pCell

    integer  :: fill( 2*int(aCell/pxReal**2), 2), rTmp( 2*int(aCell/pxReal**2), 2)
    integer  :: nFill, pxTmp
    real(b8) :: prob, r, ui, uf, w

    call getSignal( rCell, pxCell, globalSignal, localSignal)

    fill = 0
    rTmp = 0
    rTmp = rCell
    if ( a(3) == 1 ) then
        ! cell is adding a pixel
        pxTmp = pxCell + 1
        rTmp(pxTmp,1:2) = b(1:2)

        w = getWork( globalSignal, localSignal(a(4)), a(3))

        ui = getEnergy( rSim, rCell, pxCell, pCell)
        uf = getEnergy( rSim, rTmp,  pxTmp, pCell)

        prob = getProb( uf, ui, w)
    else
        ! cell is removing a pixel
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

            ui = getEnergy( rSim, rCell, pxCell, pCell)
            uf = getEnergy( rSim, rTmp,  pxTmp, pCell)

            prob = getProb( uf, ui, w)
        end if
    end if
    ! write(*,*) ' prob, w = ', prob, w
    call random_number(r)
    if ( r < prob ) then
         rCell =  rTmp
        pxCell = pxTmp
        ! write(*,*) '  success: a, b', a, b
    end if
end subroutine getElemStep


! subroutine of all the steps neccasary for one elementary initialization time-step
subroutine getItlStep( a, b, rSim, rCell, pxCell, pCell)
    implicit none
    real(b8), intent(in)    :: pCell
    integer,  intent(in)    :: a(4), b(4), rSim(2)
    integer,  intent(inout) :: rCell(:,:), pxCell
    integer  :: fill( 2*int(aCell/pxReal**2), 2), rTmp( 2*int(aCell/pxReal**2), 2)
    integer  :: nFill, pxTmp
    real(b8) :: prob, r, ui, uf, w

    fill = 0
    rTmp = 0
    rTmp = rCell
    if ( a(3) == 1 ) then
        ! cell is adding a pixel
        pxTmp = pxCell + 1
        rTmp(pxTmp,1:2) = b(1:2)

        w  = 0.0_b8
        ui = getEnergy( rSim, rCell, pxCell, pCell)
        uf = getEnergy( rSim, rTmp,  pxTmp, pCell)
        prob = getProb( uf, ui, w)
    else
        ! cell is removing a pixel
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
            w  = 0.0_b8
            ui = getEnergy( rSim, rCell, pxCell, pCell)
            uf = getEnergy( rSim, rTmp,  pxTmp, pCell)
            prob = getProb( uf, ui, w)
        end if
    end if
    call random_number(r)
    if ( r < prob ) then
         rCell =  rTmp
        pxCell = pxTmp
        ! if ( a(3) == 1 ) then
        !     write(*,*) ' cell +1 :', b(1:2)
        ! else
        !     write(*,*) ' cell -1 :', b(1:2)
        ! end if
    end if
end subroutine getItlStep


! count number of pixels neighboring ECM
subroutine perimCheck( rCell, pxCell, rSim, perim)
    ! lCell = lattice site from cell list
    implicit none
    integer, intent(in)  :: rCell(:,:), pxCell, rSim(2)
    integer, intent(out) :: perim
    integer :: i, j, k, ecmCheck, nn(2)

    perim = 0
    do i = 1, pxCell
        ! check neighbors of all cell pixels
        do j = 1, 4
            call nnGet( j, rCell(i,1:2), rSim, nn)
            ecmCheck = 1
            do k = 1, pxCell
                if ( nn(1) == rCell(k,1) .AND. nn(2) == rCell(k,2) ) then
                    ecmCheck = 0
                    exit
                end if
            enddo
            perim = perim + ecmCheck
        enddo
    enddo

end subroutine perimCheck


end module
