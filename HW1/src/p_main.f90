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

    call s_read_user_input()

    call s_initialize_global_parameters()

    call s_print_user_input()

    allocate(phi1(0:N))
    allocate(phi2(0:N))
    allocate(gamma(0:N-1))
    allocate(b(0:N))
    allocate(Q(0:N))
    allocate(A(0:N, 0:2))

    if (gamma1 == 0d0) then
        ! Linear equation
        call cpu_time(tstart)

        call s_compute_gamma(gamma, phi1)
        call s_compute_q(Q)
        call s_compute_b(b, Q)
        call s_compute_a(A, gamma)
        call s_TDMA(A, phi1, b)

        call cpu_time(tstop)
    else
        ! Nonlinear equation
        call cpu_time(tstart)

        do
            call s_compute_gamma(gamma, phi1)
            call s_compute_q(Q)
            call s_compute_b(b, Q)
            call s_compute_a(A, gamma)
            call s_TDMA(A, phi2, b)

            call s_compute_error(phi1, phi2, error)

            if (error < tol) exit

            phi1 = phi2

            n_iter = n_iter + 1
        end do

        call cpu_time(tstop)
    end if

    print*, "Execution time: ", tstop  - tstart, " s"
    print*, "Iterations: ", n_iter

    call s_write_data(phi1)

end program p_main
