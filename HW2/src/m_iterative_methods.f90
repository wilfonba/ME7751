module m_iterative_methods

    ! Dependencies ==========================================
    use m_global_parameters

    use m_helpers
    ! =======================================================

    implicit none

    private; public :: s_jacobi_iteration, &
        s_gauss_seidel_iteration, &
        s_SOR_iteration, &
        s_geometric_multigrid, &
        s_estimate_error

    integer :: i, j, k ! Generic indexe

contains

    subroutine s_jacobi_iteration()

        integer :: s, d

        s = idx(2)
        d = idx(0)

        ! Lower boundary
        do i = 1, N-1
            T(i,-1,s) = T(i,1,s)
            T(i,N+1,s) = T(i,N-1,s)
        end do

        ! Interior Points
        do i = 1, N-1
            do j = 0, N
                T(i,j,d) = h**2d0*Q(i,j)/4d0 &
                    + (1d0/4d0)*(T(i-1,j,s) + T(i+1,j,s) + T(i,j-1,s) + T(i,j+1,s))
            end do
        end do

    end subroutine s_jacobi_iteration

    subroutine s_gauss_seidel_iteration()

    end subroutine s_gauss_seidel_iteration

    subroutine s_SOR_iteration()

    end subroutine s_SOR_iteration

    subroutine s_geometric_multigrid()

    end subroutine s_geometric_multigrid

    subroutine s_estimate_error(n_iter)

        integer :: n_iter

        do i = 1, N-11
            do j = 1, N-1
                res(i,j) = abs(T(i,j,idx(0)) - T(i,j,idx(2)))
            end do
        end do

        resG(n_iter) = maxval(res)

    end subroutine s_estimate_error

end module m_iterative_methods
