module m_linear_algebra

    ! Dependencies ==========================================
    use m_global_parameters

    use m_helpers
    ! =======================================================

    implicit none

    private; public :: s_TDMA

    contains

        ! This subroutine solves the system Ax = b for a tridiagonal
        ! matrix A. Its inputs are
        !   A - an N+1 x N+1 matrix with a discretzation
        !   x - an N+1 x 1 vector to store the solution in
        !   b - an N+1 x 1 vector storing the right hand side
        subroutine s_TDMA(A, x, b)

            real(kind(0d0)), dimension(0:N, 0:2) :: A
            real(kind(0d0)), dimension(0:N) :: x, b, y
            integer :: i
            real(kind(0d0)) :: w

            ! In place LU factorization
            do i = 1, N
                w = A(i,0)/A(i-1,1)

                A(i,0) = A(i,0) - w*A(i-1,1)
                A(i,1) = A(i,1) - w*A(i-1,2)
                b(i) = b(i) - w*b(i - 1)

            end do

            ! Backward substitution
            x(N) = b(N)/A(N,1)
            do i = N-1, 0, -1
               x(i) = (b(i) - A(i,2)*x(i + 1))/A(i,1)
            end do

        end subroutine s_TDMA

end module m_linear_algebra
