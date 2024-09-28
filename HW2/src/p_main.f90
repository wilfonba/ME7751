program p_main

    ! Dependencies ==========================================
    use m_iterative_methods

    use m_global_parameters

    use m_helpers
    ! =======================================================

    implicit none

    integer :: i

    real(kind(0d0)) :: resR
    real(kind(0d0)) :: tstart, tstop
    real :: t_run

    resR = 1

    call s_read_user_input()

    call s_initialize_global_parameters()

    call s_print_user_input()

    do
        select case(solver)
            case (0) ! Jacobi
                call s_jacobi_iteration()
            case (1) ! Gauss-Seidel
                call s_gauss_seidel_iteration()
            case (2) ! SOR
                call s_SOR_iteration()
            case (3) ! Geometric Multigrid
                call s_geometric_multigrid()
        end select

        call s_estimate_error(n_iter)

        resR = resG(n_iter)/resG(1)

        if (n_iter > 0) then
            print("(A11, I4, A11, ES14.7)"), "Iteration: ", n_iter, " Residual: ", resR
        end if

        n_iter = n_iter + 1

        if (resR < tol) then
            call s_write_data(T(0:N, 0:N, idx(0)))
            print*, "Iterative method converged"
            stop
        end if

        if (n_iter == max_iter) then
            call s_write_data(T(0:N,0:N,idx(0)))
            print*, "Iterative method failed to converge"
            stop
        end if

        call s_shift_left()

    end do

end program p_main
