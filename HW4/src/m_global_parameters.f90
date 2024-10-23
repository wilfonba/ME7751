module m_global_parameters

    ! Dependencies ==========================================

    ! =======================================================

    implicit none

    real(kind(0d0)), parameter :: pi = 3.141592653589793

    ! User Inputs
    real(kind(0d0)) :: dt    ! Delta t
    integer :: Re            ! Reynolds number
    logical :: bench         ! whether to benchmark or not
    real(kind(0d0)) :: t_start ! start time
    real(kind(0d0)) :: t_stop  ! stop time
    real(kind(0d0)) :: t_save  ! save interval

    integer :: n_fields = 2

    real(kind(0d0)), dimension(:), allocatable :: xs, ys ! global x cooredinates of FD points

    public

contains

    subroutine s_initialize_global_parameters(N, h)

        integer :: N ! cells in x-direction
        integer :: i ! standard loop iterator
        real(kind(0d0)) :: h

        h = 1d0/N

        allocate(xs(0:N))
        allocate(ys(0:N))

        do i = 0,N
            xs(i) = i*h
            ys(i) = i*h
        end do

    end subroutine s_initialize_global_parameters

end module m_global_parameters
