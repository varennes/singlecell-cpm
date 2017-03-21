module probability
! module for calculating probability
use parameters
use sysconfig

contains

real(b8) function getWork( global, local, sigma)
    implicit none
    real(b8), intent(in) :: global, local
    integer,  intent(in) :: sigma

    getWork = w0 * (local-global) * (-1.0_b8)**(1-sigma)
end function getWork


! evaluate energy of the configuration
real(b8) function getEnergy( rSim, rCell, pxCell)
    implicit none
    integer,  intent(in) :: rSim(2), rCell(:,:), pxCell
    integer :: nn(2)
    integer :: i, j
    real(b8) :: areaCost, contact, sum1

    ! sum energy contribution due to J
    sum1 = 0.0_b8
    do i = 1, pxCell
        do j = 1, 4
            call nnGet( j, rCell(i,1:2), rSim, nn)
            contact = jCheck( rCell(i,1:2), nn, rCell, pxCell)
            sum1 = sum1 + contact
        enddo
    enddo
    ! energy contribution due to area
    areaCost = lArea * ( real(pxCell) - (aCell/pxReal**2))**2
    getEnergy = sum1 + areaCost

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

end module
