! Author: Lars Frogner
! approved to be borrowed by Kristine Onsum Moseid
module tridiag_module
implicit none

private
public solve_tridiag

integer, parameter :: DP = kind(0.0d0) ! Kind needed for double precision

contains

subroutine solve_tridiag(a, b, c, v)

    ! This subroutine solves a general tridiagonal system of equations.
    ! The terms in the substitutions can be confusing since some arrays
    ! are used for multiple quantities in order to preserve memory. The
    ! terms are therefore annotated with their mathematical meaning as
    ! described in the report for FYS4150 project 1.

    real(DP), intent(in)    :: a(:) ! Elements of lower diagonal
    real(DP), intent(in)    :: b(:) ! Elements of main diagonal
    real(DP), intent(in)    :: c(:) ! Elements of upper diagonal
    real(DP), intent(inout) :: v(:) ! In: elements of right hand side vector
                                    ! Out: elements of solution vector

    real(DP), allocatable :: b_tilde(:) ! For holding main diagonal elements after forward subst.
    real(DP)              :: abfrac     ! For holding reused values
    integer               :: n          ! Number of equations
    integer               :: i          ! Loop counter
    integer               :: im1        ! One step behind loop counter

    n = size(v)

    ! Make sure input arrays have the correct sizes
    if (size(a) /= n-1 .or. size(b) /= n .or. size(c) /= n-1) then
        stop 'solve_tridiag: a and c must have length n-1, and b must have length n.'
    end if

    allocate(b_tilde(n))

    b_tilde = b

    ! Forward substitution
    do i = 2, n

        im1 = i - 1

        abfrac = a(im1)/b_tilde(im1)
        !          |      |
        !       a_{i-1}  b_tilde_{i-1}

        b_tilde(i) = b(i) - c(im1)*abfrac
        !   |        |        |
        !   |       b_i      c_{i-1}
        !   b_tilde_i

        v(i) = v(i) - v(im1)*abfrac
        ! |     |        |
        ! |    s_i      s_tilde_{i-1}
        ! s_tilde_i

    end do

    ! -- Backward substitution:

    v(n) = v(n)/b_tilde(n)
    ! |     |      |
    ! v_n   |     b_tilde_n
    !    s_tilde_n

    do i = (n-1), 1, -1

        v(i) = (v(i) - c(i)*v(i+1))/b_tilde(i)
        ! |      |      |      |       |
        ! v_i    |     c_i  v_{i+1}   b_tilde_{i-1}
        !     s_tilde_i

    end do

    deallocate(b_tilde)

end subroutine solve_tridiag

end module tridiag_module
