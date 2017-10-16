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


real(b8) function getWorkAdpt1( global, local, sigma)
    implicit none
    real(b8), intent(in) :: global, local
    integer,  intent(in) :: sigma

    getWorkAdpt1 = w0 * (local-global) * (-1.0_b8)**(1-sigma) / global
end function getWorkAdpt1


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
    contact   = alpha  * real(perim)
    perimCost = lPerim * ( real(perim)  - (pCell/pxReal))**2
    ! energy contribution due to area
    areaCost  = lArea  * ( real(pxCell) - (aCell/pxReal**2))**2

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


subroutine getVectorUpdate( plrVec, rCell, pxCell, com, deltaCOM, globalSignal, localSignal, rSim)
    implicit none
    real(b8), intent(inout) :: plrVec(2), globalSignal, localSignal(:)
    real(b8), intent(in)    :: com(2), deltaCOM(2)
    integer,  intent(in)    :: rCell(:,:), pxCell, rSim(2)
    real(b8) :: conc, norm, q(2)
    integer  :: i, j, k, edgeCheck, edgeTotal, nn(2)

    call getSignal( rCell, pxCell, globalSignal, localSignal)

    q = 0.0_b8
    edgeTotal = 0
    ! identify edge pixels
    do i = 1, pxCell
        edgeCheck = 4
        ! check neighbors of all cell pixels
        do j = 1, 4
            call nnGet( j, rCell(i,1:2), rSim, nn)
            do k = 1, pxCell
                if ( nn(1) == rCell(k,1) .AND. nn(2) == rCell(k,2) ) then
                    edgeCheck = edgeCheck - 1
                end if
            enddo
        enddo
        if ( edgeCheck > 0 ) then
            edgeTotal = edgeTotal + 1
            conc = (localSignal(i) - globalSignal) / globalSignal
            norm = dsqrt( dot_product( real(rCell(i,:))-com, real(rCell(i,:))-com))
            do j = 1, 2
                q(j) = q(j) + conc * (real(rCell(i,j))-com(j)) / norm
            enddo
        end if
    enddo
    ! write(*,*) 'q = ', q, 'edgeTotal =', edgeTotal, ' qnorm =', q / real(edgeTotal)
    q = q / real(edgeTotal)
    plrVec = (1.0_b8 - rVec)*plrVec + eVec*q + nVec*deltaCOM

end subroutine getVectorUpdate


! subroutine of all the steps neccasary for one elementary time-step
! Cell has fixed polarity vector which will determine the work
subroutine getVectorStep( a, b, rSim, rCell, pxCell, pCell ,plrVec)
    implicit none
    integer,  intent(in)    :: a(4), b(4), rSim(2)
    integer,  intent(inout) :: rCell(:,:), pxCell
    real(b8), intent(in)    :: plrVec(2)
    real(b8), intent(inout) :: pCell
    integer  :: fill( 4*int(aCell/pxReal**2), 2), rTmp( 4*int(aCell/pxReal**2), 2)
    integer  :: nFill, pxTmp
    real(b8) :: prob, r, ui, uf, w
    real(b8) :: dxTmp, comNew(2), comOld(2), vec(2)

    fill = 0
    rTmp = rCell
    call getCOM( rCell, comOld)
    if ( a(3) == 1 ) then
        ! cell is adding a pixel
        pxTmp = pxCell + 1
        if ( pxTmp > 4*int(aCell/pxReal**2) ) then
            ! cell is adding more pixels than allocated
            prob = 0.0_b8
        else
            rTmp(pxTmp,1:2) = b(1:2)
            call getCOM( rTmp, comNew)
            w = dot_product( comNew-comOld, plrVec)
            dxTmp = dsqrt( dot_product( comNew-comOld, comNew-comOld) )
            if ( w /= 0.0 .AND. dxTmp > 1e-10 ) then
                w = w / dxTmp
            end if
            ui   = getEnergy( rSim, rCell, pxCell, pCell)
            uf   = getEnergy( rSim, rTmp,  pxTmp, pCell)
            prob = getProb( uf, ui, w)
        end if
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
            call getCOM( rTmp, comNew)
            w = dot_product( comNew-comOld, plrVec)
            dxTmp = dsqrt( dot_product( comNew-comOld, comNew-comOld) )
            if ( w /= 0.0 .AND. dxTmp > 1e-10 ) then
                w = w / dxTmp
            end if

            ui   = getEnergy( rSim, rCell, pxCell, pCell)
            uf   = getEnergy( rSim, rTmp,  pxTmp, pCell)
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
end subroutine getVectorStep


! subroutine of all the steps neccasary for one elementary time-step
subroutine getElemStep( a, b, rSim, rCell, pxCell, pCell, globalSignal, localSignal)
    implicit none
    integer, intent(in)     :: a(4), b(4), rSim(2)
    integer, intent(inout)  :: rCell(:,:), pxCell
    real(b8), intent(inout) :: globalSignal, localSignal(:), pCell
    integer  :: fill( 4*int(aCell/pxReal**2), 2), rTmp( 4*int(aCell/pxReal**2), 2)
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
        ! w = getWorkAdpt1( globalSignal, localSignal(a(4)), a(3))

        ui   = getEnergy( rSim, rCell, pxCell, pCell)
        uf   = getEnergy( rSim, rTmp,  pxTmp, pCell)
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
            ! w = getWorkAdpt1( globalSignal, localSignal(b(4)), a(3))

            ui   = getEnergy( rSim, rCell, pxCell, pCell)
            uf   = getEnergy( rSim, rTmp,  pxTmp, pCell)
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
    integer  :: fill( 4*int(aCell/pxReal**2), 2), rTmp( 4*int(aCell/pxReal**2), 2)
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
