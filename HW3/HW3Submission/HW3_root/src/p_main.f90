program p_main

    ! Dependencies ==========================================
    use m_linear_algebra

    use m_global_parameters

    use m_helpers
    ! =======================================================

    implicit none

    integer :: N          ! global number of cells in x-direction
    integer :: n_test = 1 ! number of tests to run
    real(kind(0d0)) :: h  ! Global grid spacing
    integer :: i, j       ! Standard loop iterator

    real(kind(0d0)), dimension(:), allocatable :: phi
    real(kind(0d0)), dimension(:), allocatable :: RHS
    real(kind(0d0)), dimension(:,:), allocatable :: A

    real(kind(0d0)) :: tr_start, tr_stop, tr_total

    integer :: n_stop  ! total # of time steps
    integer :: n_save ! time step interval at which to save

    call s_read_user_input(N)

    call s_initialize_global_parameters(N, h)

    call s_print_user_input(N)

    call s_intitialize_problem()

    n_stop = (t_stop-t_start)/dt
    n_save = t_save/dt

    if (bench) then
        n_test = 1e2
    end if

    tr_total = 0d0

    do i = 1, n_test

        call s_compute_initial_condition(phi, N)

        call s_save_data(phi, N)

        call cpu_time(tr_start)

        do j = 1, n_stop

            if (time_stepper == 0) then
                call s_compute_RHS(phi, RHS, N, h)
                phi = RHS
            else
                call s_compute_RHS(phi, RHS, N, h)
                call s_compute_A(phi, A, N, h)
                call s_tdma(A, phi, RHS, N)
            end if

            ! Saving and terminating
            if (mod(j, n_save) == 0) then
                call s_save_data(Phi, N)
                if (j == n_stop) continue
            end if

        end do

        call cpu_time(tr_stop)

        tr_total = tr_total + tr_stop - tr_start

        if (bench) then
            print*, "Test ", i, " completed."
        end if

    end do

    call s_finalize_problem()

    print*, "Runtime Inforation:"
    print*, "Time To Solution: ", tr_total/n_test
    print*, "Time/Step: ", tr_total/(n_test*n_stop)

contains

    subroutine s_intitialize_problem()

        allocate(Phi(0:n))
        allocate(RHS(0:n))

        if (time_stepper == 1 .or. time_stepper == 2) then
            allocate(A(0:n, 0:2))
        end if

        call s_open_data_file(N)

    end subroutine s_intitialize_problem

    subroutine s_finalize_problem()

        deallocate(Phi)
        deallocate(RHS)

        if (time_stepper == 1 .or. time_stepper == 3) then
            deallocate(A)
        end if

        call s_close_data_file()

    end subroutine s_finalize_problem

end program p_main
