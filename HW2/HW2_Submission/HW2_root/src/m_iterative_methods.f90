module m_iterative_methods

    ! Dependencies ==========================================
    use m_global_parameters

    use m_helpers
    ! =======================================================

    implicit none

    private; public :: s_jacobi_iteration, &
        s_gauss_seidel_iteration, &
        s_SOR_iteration, &
        s_compute_residual, &
        s_restriction, &
        s_error_solve, &
        s_prolongation, &
        s_correction

    integer :: i, j, k ! Generic indexe

contains

    ! Performs one iteration of Jacobi Iteration. Inputs are:
    ! As - Solution at iteration N
    ! Ad - Solution at iteration N + 1
    ! b - right hand side
    ! NC - number of cells in each direction
    ! hl - grid spacing
    subroutine s_jacobi_iteration(As, Ad, b, NC, hl)

        real(kind(0d0)), dimension(0:, -1:) :: As, Ad
        real(kind(0d0)), dimension(0:, 0:) :: b
        integer :: NC
        real(kind(0d0)) :: hl

        ! Lower boundary
        do i = 1, NC-1
            As(i,-1) = As(i,1)
            As(i,N+1) = As(i,N-1)
        end do

        ! Interior Points
        do i = 1, NC-1
            do j = 0, NC
                Ad(i,j) = hl**2d0*b(i,j)/4d0 &
                    + (1d0/4d0)*(As(i-1,j) + As(i+1,j) + As(i,j-1) + As(i,j+1))
            end do
        end do

    end subroutine s_jacobi_iteration

    ! Performs one iteration of Gauss-Seidel Iteration. Inputs are:
    ! As - Solution at iteration N
    ! Ad - Solution at iteration N + 1
    ! b - right hand side
    ! NC - number of cells in each direction
    ! hl - grid spacing
    subroutine s_gauss_seidel_iteration(As, Ad, b, NC, hl)

        real(kind(0d0)), dimension(0:, -1:) :: As, Ad
        real(kind(0d0)), dimension(0:, 0:) :: b
        integer :: NC
        real(kind(0d0)) :: hl

        ! Lower boundary
        do i = 1, NC-1
            Ad(i,-1) = As(i,1)
        end do

        ! Interior Points
        do i = 1, NC-1
            do j = 0, NC-1
                Ad(i,j) = hl**2d0*b(i,j)/4d0 &
                    + (1d0/4d0)*(Ad(i-1,j) + As(i+1,j) + Ad(i,j-1) + As(i,j+1))
            end do
        end do

        ! Upper boundary
        do i = 1, NC-1
            As(i,NC+1) = As(i,NC-1)
            Ad(i,NC) = hl**2d0*b(i,j)/4d0 &
                    + (1d0/4d0)*(Ad(i-1,j) + As(i+1,j) + Ad(i,j-1) + As(i,j+1))
        end do

    end subroutine s_gauss_seidel_iteration

    ! Performs one iteration of Gauss-Seidel Iteration. Inputs are:
    ! As - Solution at iteration N
    ! Ad - Solution at iteration N + 1
    ! b - right hand side
    ! NC - number of cells in each direction
    ! hl - grid spacing
    subroutine s_SOR_iteration(As, Ad, b, NC, hl)

        real(kind(0d0)), dimension(0:, -1:) :: As, Ad
        real(kind(0d0)), dimension(0:, 0:) :: b
        integer :: NC
        real(kind(0d0)) :: hl

        ! Lower boundary
        do i = 1, NC-1
            Ad(i,-1) = As(i,1)
        end do

        ! Interior Points
        do i = 1, NC-1
            do j = 0, NC-1
                 Ad(i,j) = (1d0-omega)*As(i,j) + omega*hl**2d0*b(i,j)/4d0 &
                    + (omega/4d0)*(Ad(i-1,j) + As(i+1,j) + Ad(i,j-1) + As(i,j+1))
            end do
        end do

        ! Upper boundary
        do i = 1, NC-1
            As(i,NC+1) = As(i,NC-1)
            Ad(i,NC) = (1d0-omega)*As(i,j) + omega*hl**2d0*b(i,j)/4d0 &
                    + (omega/4d0)*(Ad(i-1,j) + As(i+1,j) + Ad(i,j-1) + As(i,j+1))
        end do

    end subroutine s_SOR_iteration

    ! Computes the residual. Inputs are:
    ! D_local - Solution the residual is computed on
    ! b - right hand side
    ! NC - numer of cells in each direction
    ! hl - grid spacing
    ! resMax - largest residual
    subroutine s_compute_residual(D_local, b, NC, hl, resMax)

        real(kind(0d0)), dimension(0:,-1:) :: D_local
        real(kind(0d0)), dimension(0:,0:) :: b
        integer :: NC
        real(kind(0d0)) :: hl, resMax

        do i = 1, NC-1
            do j = 0, NC
                res(i,j) = (D_local(i+1,j) + D_local(i-1,j) - 4d0*D_local(i,j) &
                    + D_local(i,j+1) + D_local(i,j-1))/(hl**2d0) + b(i,j)
            end do
        end do

        resMax = maxval(abs(res(0:NC, 0:NC)))

    end subroutine s_compute_residual

    ! Performs a restriction from a fine to coarse mesh. Inputs are:
    ! Tf - the fine grid solution
    ! Tc - the coarse grain solution
    subroutine  s_restriction(Tf, Tc)

        real(kind(0d0)), dimension(0:,0:) :: Tf, Tc

        integer :: i, j

        ! Left and Right boundaries
        do i = 0, N/2
            Tc(0,i) = 0d0
            Tc(N/2-1,i) = 0d0
        end do

        ! Interior points
        do i = 0, N/2
            do j = 0, N/2
                Tc(i,j) = Tf(2*i, 2*j)
            end do
        end do

    end subroutine s_restriction

    ! Performs the sub-grid error solve. Inputs are:
    ! Tc - solution on coarse grid
    ! hl - grid spacing on coarse grid
    ! NC - number of cells in each direction
    subroutine s_error_solve(Tc, hl, NC)

        real(kind(0d0)), dimension(0:, 0:) :: Tc
        real(kind(0d0)) :: hl, resmax
        integer :: NC, n_iterl
        real(kind(0d0)), dimension(0:max_iter) :: resL
        real(kind(0d0)), dimension(0:NC, 0:NC) :: b

        n_iterl = 0

        res(NC,:) = 0d0
        resMax = 1d0

        b = Tc
        do
            select case(solver)
                case (0) ! Jacobiq
                    print*, "Jacobi iteration is not supported for algebraic multigrid"
                    exit
                case (1) ! Gauss-Seidel
                    call s_gauss_seidel_iteration(Tc(0:NC,-1:NC+1), Tc(0:NC,-1:NC+1), b, NC, 2*h)
                case (2) ! SOR
                    call s_SOR_iteration(Tc(0:NC,-1:NC+1), Tc(0:NC,-1:NC+1), b, NC, 2*h)
            end select

            call s_compute_residual(Tc(0:NC,-1:NC+1), b, NC, 2*h, resMax)

            if (resMax< tol) then
                exit
            end if

            n_iterl = n_iterl + 1

            if (n_iterl > max_iter) then
                print*, "Error solve failed to converge"
                stop
            end if

        end do

    end subroutine s_error_solve

    ! Perform prolongation. Inputs are:
    ! Tc - error solution on fine grid
    ! Tf - error solution on fine grid
    subroutine s_prolongation(Tc, Tf)

        real(kind(0d0)), dimension(0:N/2, 0:N/2) :: Tc
        real(kind(0d0)), dimension(0:N, 0:N) :: Tf

        integer :: i, j

        Tf = 0d0

        ! Even indices
        do i = 1, N/2-1
            do j = 0, N/2
                Tf(2*i,2*j) = Tc(i, j)
            end do
        end do

        ! x-direction
        do i = 1, N/2
            do j = 0, N/2
                Tf(2*i - 1, 2*j) = (Tf(2*(i-1),2*j) + Tf(2*i, 2*j))/2d0
            end do
        end do

        ! y-direction
        do i = 1, N/2
            do j = 1, N/2
                Tf(2*i, 2*j - 1) = (Tf(2*i,2*(j-1)) + Tf(2*i, 2*j))/2d0
            end do
        end do

        ! Centers
        do i = 1, N/2
            do j = 1, N/2
                Tf(2*i - 1, 2*j - 1) = (Tf(2*(i-1), 2*j) + Tf(2*i, 2*j) &
                        + Tf(2*(i-1),2*(j-1)) + Tf(2*i, 2*(j-1)))/4d0
            end do
        end do

    end subroutine s_prolongation

    ! Apply the correction to the fine grid solution. Inputs are:
    ! E - calculated error from multigrid step
    ! T - solution on the fine grid
    subroutine s_correction(E, T)

        real(kind(0d0)), dimension(0:N, 0:N) :: E, T
        integer :: i, j

        do i = 0,N
            do j = 0,N
                T(i,j) = T(i,j) + E(i,j)
            end do
        end do

    end subroutine s_correction

end module m_iterative_methods

