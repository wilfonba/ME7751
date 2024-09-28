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
    integer, dimension(0:2) :: idx ! indicies of past steps
    real(kind(0d0)), dimension(0:2) :: z ! z for error estimationg

    public

contains

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

        ! Apply left and right boundary conditions
        do i = 0, N
            T(0,i,0) = 2d0*ys(i)**3d0 - 3d0*ys(i)**2d0 + 1
            T(0,i,1) = T(0,i,0)
            T(0,i,2) = T(0,i,0)

            T(N,i,0) = 0d0
            T(N,i,1) = 0d0
            T(N,i,2) = 0d0
        end do

        ! Apply initial guess
        do i = 1, N-1
            do j = 0, N-1
                T(i,j,0) = 0d0
            end do
        end do

    end subroutine s_initialize_global_parameters

end module m_global_parameters
