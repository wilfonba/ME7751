module m_helpers

    ! Dependencies ==========================================
    use m_global_parameters
    ! =======================================================

    implicit none

    private; public :: s_compute_initial_condition, &
        s_read_user_input, &
        s_print_user_input, &
        s_compute_RHS, &
        s_compute_A, &
        s_print_2d_array, &
        s_open_data_file, &
        s_save_data, &
        s_close_Data_file

contains

    ! This subroutine assigns the initial condition for the
    ! nonlinear solve. It's arguments are:
    !   Phi - an N + 1 x 1 array storing phi
    !   N - the length of Phi
    subroutine s_compute_initial_condition(Phi, N)

        real(kind(0d0)), dimension(0:) :: Phi
        integer :: i, N

        do i = 0, N
            Phi(i) = (4d-1*pi)**(-5d-1)*exp(-2.5d0*(xs(i) - 10)**2d0)
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

        namelist /user_inputs/ gamma, U, bench, time_stepper, N, dt, &
                                t_start, t_stop, t_save

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

        print*, "Solving the 1D Steady Transport Equation with:"
        print("(A14, I1)"), "time_stepper: ", time_stepper
        print("(A14, F5.3)"), "gamma: ", gamma
        print("(A14, F5.3)"), "U: ", U
        print("(A14, I5)"), "N: ", N
        print("(A14, L1)"), "bench: ", bench
        print("(A14, ES10.4)"), "dt: ", dt

    end subroutine s_print_user_input

    subroutine s_compute_RHS(phi, RHS, N, h)

        real(kind(0d0)), dimension(0:) :: phi
        real(kind(0d0)), dimension(0:) :: RHS
        integer :: N
        integer :: i ! < default iteratora
        real(kind(0d0)) :: h, c, d

        d = gamma*dt/(h**2d0)
        c = U*dt/h

        if (time_stepper == 0) then ! Explicit Euler

            ! Left BC
            RHS(0) = 0d0

            ! Interior Points
            do i = 1, N-1
                RHS(i) = (d + c/2)*phi(i - 1) + (1 - 2d0*d)*phi(i) + (d - c/2d0)*phi(i + 1)
            end do

            ! Right BC
            RHS(N) = 0d0

        elseif (time_Stepper == 1) then ! Implicit Euler

            ! Left BC
            RHS(0) = 0d0

            ! Interior Points
            do i = 1, N-1
                RHS(i) = phi(i)
            end do

            ! Right BC
            RHS(N) = 0d0

        elseif (time_stepper == 2) then ! Crank.Nicolson

           ! Left BC
            RHS(0) = 0d0

            ! Interior Points
            do i = 1, N-1
                RHS(i) = (d/2d0 + c/4d0)*phi(i - 1) + (1 - d)*phi(i) + &
                         (d/2d0 - c/4d0)*phi(i + 1)
            end do

            ! Right BC
            RHS(N) = 0d0

        else
            print*, "Invalid time_stepper. Exiting..."
            call exit(1)
        end if

    end subroutine s_compute_RHS

    subroutine s_compute_A(phi, A, N, h)

        real(kind(0d0)), dimension(0:) :: phi
        real(kind(0d0)), dimension(0:, 0:) :: A
        integer :: N
        integer :: i ! generic iterator
        real(kind(0d0)) :: h, c, d

        d = gamma*dt/(h**2d0)
        c = U*dt/h

        if (time_stepper == 0) then ! Expliti Euler
            ! No A required
        elseif (time_Stepper == 1) then ! Implicit Euler

            ! First row
            A(0,0:2) = [0, 1, 0]

            ! Interior rows
            do i = 1, N-1
                A(i,:) = [-c/2d0 - d, 1 + 2d0*d, c/2d0 - d]
            end do

            ! Last row
            A(N,0:2) = [0, 1, 0]

        elseif (time_stepper == 2) then ! Crank-Nicolson

            ! First row
            A(0,0:2) = [0, 1, 0]

            ! Interior rows
            do i = 1, N-1
                A(i,:) = [-c/4d0 - d/2d0, 1 + d, c/4d0 - d/2d0]
            end do

            ! Last row
            A(N,0:2) = [0, 1, 0]

        else
            print*, "Invalid time_stepper. Exiting..."
            call exit(1)
        end if

    end subroutine s_compute_A

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

    subroutine s_open_data_file(N)

        integer :: N, i

        character(LEN=100) :: file_name = 'output.csv'

        open (3, FILE=trim(file_name), &
              FORM='formatted', &
              STATUS='replace')

        do i = 0, N-1
            write(3,"(F16.7, A2)",advance="no") xs(i), ","
        end do
        write(3, "(F16.7)") xs(N)

    end subroutine s_open_data_file

    subroutine s_save_data(phi, N)

        real(kind(0d0)), dimension(0:) :: phi
        integer :: i, N

        do i = 0, N-1
            write(3,"(F16.7, A2)",advance="no") phi(i), ","
        end do
        write(3, "(F16.7)") phi(N)

    end subroutine s_save_data

    subroutine s_close_data_file()

        close(3)

    end subroutine s_close_data_file

end module m_helpers

