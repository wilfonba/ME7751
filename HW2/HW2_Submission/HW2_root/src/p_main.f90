program p_main

    ! Dependencies ==========================================
    use m_iterative_methods

    use m_global_parameters

    use m_helpers
    ! =======================================================

    implicit none

    integer :: i, j

    real(kind(0d0)) :: resR
    real(kind(0d0)) :: tstart, tstop, ttotal
    integer :: n_run = 1;
    real(kind(0d0)) :: resMax

    call s_read_user_input()

    if (bench) then
        n_run = 1e4/N
    end if

    call s_initialize_global_parameters()

    call s_print_user_input()

    if (.not. multigrid) then

        do i = 1, n_run

            call s_initial_guess()

            n_iter = 0
            resR = 1
            resG(1) = 1

            if (bench) then
                call cpu_time(tstart)
            end if

            do
                select case(solver)
                    case (0) ! Jacobi
                        call s_jacobi_iteration(T(:,:,idx(2)), T(:,:,idx(0)), Q, N, h)
                    case (1) ! Gauss-Seidel
                        call s_gauss_seidel_iteration(T(:,:,idx(2)), T(:,:,idx(0)), Q, N, h)
                    case (2) ! SOR
                        call s_SOR_iteration(T(:,:,idx(2)), T(:,:,idx(0)), Q, N, h)
                end select

                call s_compute_residual(T(:,:,idx(0)), Q, N, h, resMax)

                resG(n_iter) = resMax

                if (n_iter > 0 .and. .not. bench) then
                    print("(A11, I6, A11, ES14.7)"), "Iteration: ", n_iter, " Residual: ", resG(n_iter)
                end if

                if (resG(n_iter) < tol) then
                    if (bench) then
                        call cpu_time(tstop)
                        ttotal = ttotal + tstop - tstart
                    end if

                    call s_write_data(T(0:N, 0:N, idx(0)))

                    print*, "Iterative method converged"
                    exit
                end if

                n_iter = n_iter + 1

                if (n_iter == max_iter) then
                    call s_write_data(T(0:N,0:N,idx(0)))
                    print*, "Iterative method failed to converge"
                    stop
                end if

                call s_shift_left()

            end do

        end do

    else

        do i=  1, n_run

            call s_initial_guess()

            n_iter = 0
            resR = 1
            resG(1) = 1

            do
                call s_apply_bcs()

                if (bench) then
                    call cpu_time(tstart)
                end if

                do j = 0, 5
                     select case(solver)
                        case (0) ! Jacobi
                            print*, "Multigrid not implemented for Jacobi iteration"
                            stop
                        case (1) ! Gauss-Seidel
                            call s_gauss_seidel_iteration(T(:,:,idx(2)), T(:,:,idx(0)), Q, N, h)
                        case (2) ! SOR
                            call s_SOR_iteration(T(:,:,idx(2)), T(:,:,idx(0)), Q, N, h)
                    end select
                    call s_shift_left()
                end do

                call s_compute_residual(T(:,:,idx(0)), Q, N, h, resMax)

                resG(n_iter) = resMax

                if (n_iter > 0 .and. .not. bench) then
                    print("(A11, I6, A11, ES14.7)"), "Iteration: ", n_iter, " Residual: ", resG(n_iter)
                end if

                if (resG(n_iter) < tol) then
                    if (bench) then
                        call cpu_time(tstop)
                        ttotal = ttotal + tstop - tstart
                    end if
                    call s_write_data(T(0:N, 0:N, idx(0)))

                    print*, "Iterative method converged"
                    exit
                end if

                n_iter = n_iter + 1

                if (n_iter == max_iter) then
                    call s_write_data(T(0:N,0:N,idx(0)))
                    print*, "Iterative method failed to converge"
                    stop
                end if

                call s_restriction(res(0:N,0:N), T(0:N/2, 0:N/2, idx(1)))

                call s_error_solve(T(0:N/2, 0:N/2, idx(1)), 2*h, N/2)

                call s_prolongation(T(0:N/2, 0:N/2, idx(1)), T(0:N, 0:N, idx(2)))

                call s_correction(T(0:N, 0:N, idx(2)), T(0:N,0:n,idx(0)))

                call s_shift_left()
            end do

        end do

    endif

    if (bench) then
        print*, "N_iter: ", n_iter
        print*, "Average runtime: ", ttotal/n_run
        print*, "Average iteration time: ", ttotal/(n_run*n_iter)
    end if

    call s_finalize_global_parameters()

end program p_main
