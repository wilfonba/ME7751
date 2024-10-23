module m_helpers

    ! Dependencies ==========================================
    use m_global_parameters
    ! =======================================================

    implicit none

    private; public :: s_compute_initial_condition, &
        s_read_user_input, &
        s_print_user_input, &
        s_print_2d_array

contains

    ! This subroutine assigns the initial condition for the
    ! nonlinear solve. It's arguments are:
    !   Q - an N + 1 x N + 1 x N + 1 array storing the solution
    !   N - the length of Phi
    subroutine s_compute_initial_condition(Q, N)

        real(kind(0d0)), dimension(0:,0:,1:) :: Q
        integer :: i, j, k, N

        do k = 1, 2
            do j = 0, N
                do i = 0, N
                    Q(i,j,k) = ((-1)**k)*xs(i)*ys(j)
                end do
            end do
        end do

    end subroutine s_compute_initial_condition

    ! This subroutine reads a namelist of user inputs.
    ! It requires no arguments
    subroutine s_read_user_input(N)

        integer :: N

        character(LEN=100) :: line
        character(LEN=100) :: file_path = 'main.inp'
        integer :: iostatus
        logical :: file_exist

        namelist /user_inputs/ N, dt, Re, t_start, t_stop, t_save, bench

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
    subroutine s_print_user_input(N)

        integer :: N

        print*, "Solving the 2D Lid Driven Cavity Problem With:"
        print("(A14, I5)"), "Re: ", Re
        print("(A14, I5)"), "N: ", N
        print("(A14, L1)"), "bench: ", bench
        print("(A14, ES10.4)"), "dt: ", dt

    end subroutine s_print_user_input

    ! Prints a formatted 2D array to the terminal for debuggin purposes. Inputs
    ! are:
    ! A - the 2D array to print
    ! div - an optionial constant to divide each number by
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

