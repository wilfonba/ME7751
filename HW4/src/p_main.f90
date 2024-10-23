program p_main

    ! Dependencies ==========================================
    use m_linear_algebra

    use m_global_parameters

    use m_helpers

    use m_vtk
    ! =======================================================

    implicit none

    integer :: N          ! global number of cells in x-direction
    integer :: n_test = 1 ! number of tests to run
    real(kind(0d0)) :: h  ! Global grid spacing
    integer :: i, j       ! Standard loop iterator

    real(kind(0d0)) :: tr_start, tr_stop, tr_total

    integer :: n_stop  ! total # of time steps
    integer :: n_save ! time step interval at which to save

    real(kind(0d0)), dimension(:,:,:), allocatable :: Q

    character(len=10) :: name

    call s_read_user_input(N)

    call s_initialize_global_parameters(N, h)

    call s_print_user_input(N)

    call s_intitialize_problem()

    call s_compute_initial_condition(Q, N)

    call s_open_vtk_data_file(N, 0)
    name = "variable0"
    call s_write_variable_to_vtk_file(Q(0:N,0:N,1),N,name)
    name = "variable1"
    call s_write_variable_to_vtk_file(Q(0:N,0:N,2),N,name)
    call s_close_vtk_data_file()

    n_stop = (t_stop-t_start)/dt
    n_save = t_save/dt

    if (bench) then
        n_test = 1e2
    end if

    tr_total = 0d0

    do i = 1, n_test

        call cpu_time(tr_start)

        do j = 1, n_stop

            ! Saving and terminating
            if (mod(j, n_save) == 0) then

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

        allocate(Q(0:N,0:N,1:2))

    end subroutine s_intitialize_problem

    subroutine s_finalize_problem()


    end subroutine s_finalize_problem

end program p_main
