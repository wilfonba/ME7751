program p_main

    ! Dependencies ==========================================
    use m_linear_algebra

    use m_global_parameters

    use m_helpers

    use m_vtk

    use m_finite_difference
    ! =======================================================

    implicit none

    integer :: N          ! global number of cells in x-direction
    integer :: n_test = 1 ! number of tests to run
    real(kind(0d0)) :: h  ! Global grid spacing
    integer :: i, j, k    ! Standard loop iterator

    real(kind(0d0)) :: tr_start, tr_stop, tr_total

    integer :: max_iter, save_iter, n_iter

    integer :: save_count = 0

    ! Scalar field of solution variables and variable updates
    ! Q(1)%sf stores the density field
    ! Q(2)%sf stores the x-velocity field
    ! Q(3)%sf stores the y-velocity field
    type(scalar_field), allocatable, dimension(:) :: Q, dUP, dU, RHS

    ! Sparse matrix storage for LHS
    type(scalar_field), allocatable, dimension(:) :: LHS

    character(len=10) :: name

    call s_read_user_input(N, max_iter, save_iter)

    call s_initialize_global_parameters(N, h)

    call s_print_user_input(N)

    call s_initialize_problem()

    if (bench) then
        n_test = 1e2
    end if

    tr_total = 0d0

    do i = 1, n_test

        call s_compute_initial_condition(Q, N)

        call s_save_data(Q, N, save_count)

        n_iter = 0

        call cpu_time(tr_start)

        do j = 1, max_iter

            call s_compute_rhs(Q, RHS, N, h)

            !call s_compute_lhs(Q, LHS, 1, N, h)

            do k = 1, 3
                call s_TDMA(LHS(i)%sf, dUP(i)%sf, RHS(i)%sf, N)
            end do

            !call s_compute_lhs(Q, LHS, 2, N, h)

            do k = 1, 3
                call s_TDMA(LHS(i)%sf, dU(i)%sf, dUP(i)%sf, N)
            end do

            call s_update_solution(Q, dU, N)

            n_iter = n_iter + 1

        end do

        call s_save_data(Q, N, 1)

        call cpu_time(tr_stop)

        tr_total = tr_total + tr_stop - tr_start

        if (bench) then
            !!print*, "Test ", i, " completed."
        end if

    end do

    print*, "Runtime Inforation:"
    print*, "Time To Solution: ", tr_total/n_test
    print*, "Time/Step: ", tr_total/(n_test*j)

    call s_finalize_problem()

contains

    subroutine s_initialize_problem()

        integer :: i

        allocate(Q(1:3))
        allocate(dU(1:3))
        allocate(DUP(1:3))
        allocate(RHS(1:3))
        allocate(LHS(1:3))

        do i = 1,3
            allocate(Q(i)%sf(0:N, 0:N))
            allocate(dU(i)%sf(0:N, 0:N))
            allocate(dUP(i)%sf(0:N, 0:N))
            allocate(RHS(i)%sf(0:N, 0:N))
            allocate(LHS(i)%sf(0:N, 0:2))
        end do

    end subroutine s_initialize_problem

    subroutine s_finalize_problem()

        integer :: i

        do i = 1,3
            deallocate(Q(i)%sf)
            deallocate(dU(i)%sf)
            deallocate(dUP(i)%sf)
            deallocate(RHS(i)%sf)
            deallocate(LHS(i)%sf)
        end do

        deallocate(Q)
        deallocate(dU)
        deallocate(DUP)
        deallocate(RHS)
        deallocate(LHS)


    end subroutine s_finalize_problem

end program p_main
