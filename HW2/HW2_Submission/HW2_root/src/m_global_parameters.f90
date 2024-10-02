module m_global_parameters

    ! Dependencies ==========================================

    ! =======================================================

    implicit none

    ! User Inputs
    integer :: N ! number if discretization points
    integer :: solver ! Which iteration type to use
    real(kind(0d0)) :: tol ! Convergence tolerance
    logical :: bench ! Benchmark or not

    integer :: n_iter, max_iter
    real(kind(0d0)), allocatable, dimension(:) :: rho ! eigenvalue estimates

    real(kind(0d0)) :: h ! Grid spacing
    real(kind(0d0)), allocatable, dimension(:) :: xs, ys ! x and y values
    real(kind(0d0)), allocatable, dimension(:,:,:) :: T ! Temperature
    real(kind(0d0)), allocatable, dimension(:,:) :: res
    real(kind(0d0)), allocatable, dimension(:,:) :: Q ! Source term
    real(kind(0d0)), allocatable, dimension(:) :: resG
    real(kind(0d0)) :: omega
    integer, dimension(0:2) :: idx ! indicies of past steps
    real(kind(0d0)), dimension(0:2) :: z ! z for error estimationg
    logical :: multigrid

    public

contains

    ! Initializes variables. Requires no input
    subroutine s_initialize_global_parameters()

        integer :: i, j ! standard loop iterator

        idx = [0, 1, 2]

        h = 1d0/N

        allocate(xs(0:N))
        allocate(ys(0:N))

        do i = 0,N
            xs(i) = i*h
            ys(i) = i*h
        end do

        allocate(rho(0:max_iter))

        allocate(T(0:N, -1:N+1, 0:2))
        allocate(res(0:N, 0:N))
        allocate(resG(0:max_iter))
        allocate(Q(0:N, 0:N))

        ! Compute Q and make initial guess
        do i = 0, N
            do j = 0, N
                Q(i,j) = 4d0*ys(j)**3d0 - 6d0*ys(j)**2d0 + 2 - 6d0*(1 - xs(i)**2d0)*(2*ys(j) - 1)
            end do
        end do

    end subroutine s_initialize_global_parameters

    ! Deallocates memory and finalizes things. Requires no inputs
    subroutine s_finalize_global_parameters()

        deallocate(xs)
        deallocate(ys)
        deallocate(rho)
        deallocate(T)
        deallocate(res)
        deallocate(resG)
        deallocate(Q)

    end subroutine s_finalize_global_parameters

end module m_global_parameters
