module m_helpers

    ! Dependencies ==========================================
    use m_global_parameters
    ! =======================================================

    implicit none

    private; public :: s_read_user_input, &
        s_print_user_input, &
        s_write_data, &
        s_shift_left, &
        s_print_2d_array

contains

    ! This subroutine reads a namelist of user inputs.
    ! It requires no arguments
    subroutine s_read_user_input()

        character(LEN=100) :: line
        character(LEN=100) :: file_path = 'main.inp'
        integer :: iostatus
        logical :: file_exist

        namelist /user_inputs/ N, solver, bench, tol, max_iter

        inquire (FILE=trim(file_path), EXIST=file_exist)

        if (file_exist) then
            open (1, FILE=trim(file_path), &
                  FORM='formatted', &
                  ACTION='read', &
                  STATUS='old')
            read (1, NML=user_inputs, iostat=iostatus)

            if (iostatus /= 0) then
                backspace (1)
                read (1, fmt='(A)') line
                print *, 'Invalid line in namelist: '//trim(line)
                call abort()
            end if

            close (1)
        else
            print*, "Missing input file"
            call abort()
        end if

    end subroutine s_read_user_input

    ! This subroutine prints the user inputs to the terminal for
    ! verification
    subroutine s_print_user_input()

        print*, "Solving the 1D Steady Transport Equation with:"
        print*, "N = ", N
        print*, "Bench ", bench
        print*, "tol = ", tol

        select case (solver)
            case (0)
                print*, "Solver = Jacobi Iteration"
            case (1)
                print*, "Solver = Gauss-Seidel Iteration"
            case (2)
                print*, "Solver = SOR Iteration"
            case (3)
                print*, "Solver = Geometric Multigrid"
            case default
                print*, "Invalid solver specified. Exiting..."
                stop
        end select

    end subroutine s_print_user_input

    ! This function writes the solution and residuals to a .csv file.
    ! Its inputs are:
    !   Phi - an N x 1 solutio nvector
    subroutine s_write_data(T_save)

        real(kind(0d0)), dimension(0:N, 0:N) :: T_save
        integer :: i, j

        open(unit=1, file="T.csv")

        ! write xs
        do j = 0, N-1
            write(1,"(F8.4, A2)", advance="no") xs(j), ", "
        end do
        write(1, "(F8.4)") xs(N)

        ! write ys
        do j = 0, N-1
            write(1,"(F8.4, A2)", advance="no") xs(j), ", "
        end do
        write(1, "(F8.4)") xs(N)

        ! Write temperature
        do i = 0, N
            do j = 0, N-1
                write(1,"(F8.4, A2)", advance="no") T_save(i,j), ", "
            end do
            write(1,"(F8.4)") T_save(i,N)
        end do

        close(1)

        open(unit=1, file="residuals.csv")
            do j = 4, n_iter
                write(1,*) rho(j)/rho(4)
            end do
        close(1)

    end subroutine s_write_data

    subroutine s_shift_left

        integer :: temp

        temp = idx(0)
        idx(0) = idx(1)
        idx(1) = idx(2)
        idx(2) = temp

    end subroutine s_shift_left

   subroutine s_print_2D_array(A, div)

        real(kind(0d0)), dimension(:, :), intent(in) :: A
        real, optional, intent(in) :: div

        integer :: i, j
        integer :: m, n
        real :: c

        m = size(A, 1)
        n = size(A, 2)

        if (present(div)) then
            c = div
        else
            c = 1
        end if

        print *, m, n

        do i = 1, m
            do j = 1, n
                write (*, fmt="(F12.4)", advance="no") A(i, j)/c
            end do
            write (*, fmt="(A1)") " "
        end do
        write (*, fmt="(A1)") " "

    end subroutine 

end module m_helpers

