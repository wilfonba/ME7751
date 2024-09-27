program p_main

    ! Dependencies ==========================================
    use m_linear_algebra

    use m_global_parameters

    use m_helpers
    ! =======================================================

    implicit none

    integer :: i

    real(kind(0d0)), allocatable, dimension(:) :: phi1, phi2, gamma, b, Q
    real(kind(0d0)), allocatable, dimension(:,:) :: A

    real(kind(0d0)) :: error
    real(kind(0d0)) :: tstart, tstop
    integer :: n_iter = 0
    integer :: n_test = 1
    real :: t_run
    real(kind(0d0)), dimension(2) :: delta

    call s_read_user_input()

    call s_initialize_global_parameters()

    call s_print_user_input()

    allocate(phi1(0:N))
    allocate(phi2(0:N))
    allocate(gamma(0:N-1))
    allocate(b(0:N))
    allocate(Q(0:N))
    allocate(A(0:N, 0:2))

    if (bench) then
        n_test = 1e6
    end if

    if (gamma1 == 0d0) then
        ! Linear equation
        do i = 1, n_test
            call cpu_time(tstart)

            call s_compute_gamma(gamma, phi1)
            call s_compute_q(Q)
            call s_compute_b(b, Q)
            call s_compute_a(A, gamma)
            call s_TDMA(A, phi1, b)

            call cpu_time(tstop)

            t_run = t_run + (tstop - tstart)
        end do
    else
        ! Nonlinear equation
        do i = 1, n_test
            n_iter = 0d0
            call s_compute_initial_guess(phi1)

            error = 1

            call cpu_time(tstart)

            do
                call s_compute_gamma(gamma, phi1)
                call s_compute_q(Q)
                call s_compute_b(b, Q)
                call s_compute_a(A, gamma)
                call s_TDMA(A, phi2, b)

                call s_estimate_error(phi1, phi2, delta, n_iter, error)

                n_iter = n_iter + 1

                if (error < tol) exit

                phi1 = phi2

            end do

            call cpu_time(tstop)

            t_run = t_run + (tstop - tstart)
        end do
    end if

    ! I/O
    call s_write_data(phi1)

    ! Post comment and memory cleanup
    print*, "Execution time: ", t_run/n_test, " s"
    print*, "Iterations: ", n_iter

    deallocate(phi1)
    deallocate(phi2)
    deallocate(gamma)
    deallocate(b)
    deallocate(Q)
    deallocate(A)

end program p_main
