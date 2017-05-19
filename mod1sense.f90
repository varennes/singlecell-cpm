module sensing
! module for chemical sensing
use parameters
use sysconfig

contains

! mean chemical concentration, linear profile
real(b8) function cLin( x)
    implicit none
    real(b8), intent(in) :: x
    cLin = c0 + (g * x / pxReal)
end function cLin


! sample local concentration measurement at each pixel
! calculate global, whole cell averaged signal
subroutine getSignal( rCell, pxCell, globalSignal, localSignal)
    implicit none
    integer,  intent(in)  :: rCell(:,:), pxCell
    real(b8), intent(out) :: globalSignal, localSignal(:)
    real(b8)  :: c, x
    integer   :: i

    localSignal(:) = 0.0_b8
    do i = 1, pxCell
        x = real(rCell(i,1))
        c = cLin(x)
        ! sample local signal from a distribution
        if( c < 100.0 )then
            call poissonrand( c, localSignal(i))
        else
            localSignal(i) = normal( c, sqrt(c))
        endif

        if( localSignal(i) < 0.0 )then
            localSignal(i) = 0.0_b8
        endif
    enddo
    globalSignal = sum( localSignal) / real( pxCell)
end subroutine getSignal


subroutine updatePolarity( pVec, comCell, rSim, rCell, pxCell, globalSignal, localSignal)
    implicit none
    real(b8), intent(inout) :: pVec(2), globalSignal, localSignal(:)
    real(b8), intent(in)    :: comCell(2)
    integer, intent(in)     :: rSim(2), rCell(:,:), pxCell
    integer  :: i, j, k, edgeCheck, edgeN, nn(2)
    real(b8) :: norm, edgeVec(2), tempVec(2)

    call getSignal( rCell, pxCell, globalSignal, localSignal)

    ! iterate over all pixels and identify edge pixels
    edgeN = 0
    edgeVec(:) = 0.0_b8
    do i = 1, pxCell
        edgeCheck = 4
        ! check neighbors of all cell pixels
        do j = 1, 4
            call nnGet( j, rCell(i,1:2), rSim, nn)
            do k = 1, pxCell
                if ( nn(1) == rCell(k,1) .AND. nn(2) == rCell(k,2) ) then
                    edgeCheck = edgeCheck - 1
                    exit
                end if
            enddo
        enddo
        if ( edgeCheck > 0 ) then
            edgeN = edgeN + 1
            ! pixel is an edge pixel
            tempVec(:) = 0.0_b8
            norm  = dsqrt( (real(rCell(i,1))-comCell(1))**2 + (real(rCell(i,2))-comCell(2))**2 )
            do j = 1, 2
                tempVec(j) = ( real(rCell(i,j)) - comCell(j) ) / norm * (localSignal(i) - globalSignal) / globalSignal
                edgeVec(j) = edgeVec(j) + tempVec(j)
            enddo
        end if
    enddo
    edgeVec = edgeVec / real(edgeN)

    pVec = (1.0_b8 - vecR)*pVec + (vecR * vecE)*edgeVec

end subroutine updatePolarity


! returns random number between 0 - 1
function ran1()
    implicit none
    real(b8) ran1,x
    call random_number(x) ! built in fortran 90 random number function
    ran1=x
end function ran1


! returns a normal distribution
function normal(mean,sigma)
    implicit none
    real(b8) normal,tmp
    real(b8) mean,sigma
    integer flag
    real(b8) fac,gsave,rsq,r1,r2
    save flag,gsave
    data flag /0/
    if (flag.eq.0) then
    rsq=2.0_b8
        do while(rsq.ge.1.0_b8.or.rsq.eq.0.0_b8) ! new from for do
            r1=2.0_b8*ran1()-1.0_b8
            r2=2.0_b8*ran1()-1.0_b8
            rsq=r1*r1+r2*r2
        enddo
        fac=sqrt(-2.0_b8*log(rsq)/rsq)
        gsave=r1*fac
        tmp=r2*fac
        flag=1
    else
        tmp=gsave
        flag=0
    endif
    normal=tmp*sigma+mean
    return
end function normal


! returns poisson distributed random number
subroutine poissonrand( lambda, k)
    ! lambda is the mean
    ! k is the random number output
    implicit none
    real(b8), intent(in)  :: lambda
    real(b8), intent(out) :: k
    real(b8) :: l, p, u

    l = exp(-1.0*lambda)
    k = 0.0
    p = 1.0

    do while( p > l )
        k = k + 1.0
        call random_number(u)
        p = p * u
    enddo

    k = k - 1.0

end subroutine poissonrand


subroutine gauss_1(a,b,x,n)
    !============================================================
    ! Solutions to a system of linear equations A*x=b
    ! Method: the basic elimination (simple Gauss elimination)
    ! Alex G. November 2009
    !-----------------------------------------------------------
    ! input ...
    ! a(n,n) - array of coefficients for matrix A
    ! b(n) - vector of the right hand coefficients b
    ! n - number of equations
    ! output ...
    ! x(n) - solutions
    ! comments ...
    ! the original arrays a(n,n) and b(n) will be destroyed
    ! during the calculation
    !===========================================================
    implicit none
    integer n
    real(b8) a(n,n), b(n), x(n)
    real(b8) c
    integer i, j, k
    !step 1: forward elimination
    !    print*,'hi'
    do k=1, n-1
        do i=k+1,n
            c=a(i,k)/a(k,k)
            a(i,k) = 0.0
            b(i)=b(i)- c*b(k)
            do j=k+1,n
                a(i,j) = a(i,j)-c*a(k,j)
            end do
        end do
    end do
    !step 2: back substitution
    x(n) = b(n)/a(n,n)
    do i=n-1,1,-1
        c=0.0
        do j=i+1,n
            c= c + a(i,j)*x(j)
        end do
        x(i) = (b(i)- c)/a(i,i)
    end do
end subroutine gauss_1


end module
