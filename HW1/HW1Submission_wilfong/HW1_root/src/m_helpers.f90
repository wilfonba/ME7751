module m_helpers

    ! Dependencies ==========================================
    use m_global_parameters
    ! =======================================================

    implicit none

    private; public :: s_compute_b, &
        s_compute_gamma, &
        s_compute_initial_guess, &
        s_read_user_input, &
        s_compute_a, &
        s_print_user_input, &
        s_compute_q, &
        s_write_data, &
        s_estimate_error

    contains

        ! This subroutine populates the RHS b. It's arguments are:
        !   b - an N+1 x 1 array to store b in
        !   gamma - an N+2 x 1 array storing gamma
        !   Phi - an N+1 x 1 array storing Phi_i
        !   Q - an N+1 x 1 arrayu storing Q
        subroutine s_compute_b(b, Q)

            real(kind(0d0)), dimension(0:N) :: b, Q
            integer :: i

            ! i = 0 RHS with BC
            b(0) = 1d0

            ! i = 1, N - 1
            do i = 1,N-1
                b(i) = Q(i)
            end do

            ! i = N RHS with BC
            b(N) = 0d0

        end subroutine s_compute_b

        ! This subroutine populates the gamma vector. It's
        ! arguments are:
        !   gamma - an N+2 x 1 array storing gamma
        !   phi - an N+1 x 1 array storing phi
        subroutine s_compute_gamma(gamma, phi)

            real(kind(0d0)), dimension(0:N - 1) :: gamma
            real(kind(0d0)), dimension(0:N) :: phi
            real(kind(0d0)) :: avg
            integer :: i

            ! intermediate gammas
            do i = 0,N-1
                avg = (phi(i + 1) + phi(i))/2d0
                gamma(i) = gamma0 + gamma1*avg
            end do

        end subroutine s_compute_gamma

        ! This subroutine calculates the initial guess for the
        ! nonlinear solve. It's arguments are:
        !   Phi - an N+1 x 2 array to store Phi_i
        subroutine s_compute_initial_guess(Phi)

            real(kind(0d0)), dimension(0:N) :: Phi
            integer :: i

            do i = 0, N
                Phi(i) = 1 - xs(i)/L
            end do

        end subroutine s_compute_initial_guess

        ! This subroutine computes the matrix A. It's arguments
        ! are:
        !   A - an N+1 x 3 array storing the diagonals of A
        !   gamma - an N+1 x 1 array storing gamma
        subroutine s_compute_a(A, gamma)

            real(kind(0d0)), dimension(0:N, 0:2) :: A
            real(kind(0d0)), dimension(0:N - 1) :: gamma
            integer :: i

            ! Lower diagonal
            do i = 1,N - 1
                A(i,0) = -(U/(2*dx) + gamma(i - 1)/(dx**2d0))
            end do

            A(N,0) = 0d0

            ! Diagonal
            a(0,1) = 1d0

            do i = 1, N-1
                A(i,1) = (gamma(i - 1) + gamma(i))/(dx**2d0)
            end do

            a(N, 1) = 1d0

            ! Upper diagonal
            A(0,2) = 0d0

            do i = 1, N-1
                A(i,2) = U/(2d0*dx) - gamma(i)/(dx**2d0)
            end do

        end subroutine s_compute_a

        ! This subroutine calculates the source term Q. Its
        ! inputs are:
        !   Q - an N x 1 array storing Q
        subroutine s_compute_q(Q)

            real(kind(0d0)), dimension(0:N) :: Q
            integer :: i

            do i = 0, N
                Q(i) = Q0 + Q1*xs(i)
            end do

        end subroutine s_compute_q

        ! This subroutine reads a namelist of user inputs.
        ! It requires no arguments
        subroutine s_read_user_input()

            character(LEN=100) :: line
            character(LEN=100) :: file_path = 'main.inp'
            integer :: iostatus
            logical :: file_exist

            namelist /user_inputs/ gamma0, gamma1, N, L, Q0, Q1, U, tol, bench

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
            print*, "gamma0: ", gamma0
            print*, "gamma1: ", gamma1
            print*, "U: ", U
            print*, "Q0: ", Q0
            print*, "Q1: ", Q1
            print*, "N: ", N

        end subroutine s_print_user_input

        ! This function writes the solution to a .csv file.
        ! Its inputs are:
        !   Phi - an N x 1 solutio nvector
        subroutine s_write_data(Phi)

            real(kind(0d0)), dimension(0:N) :: Phi
            integer :: i

            open(unit=1, file="output.csv")

            do i = 0, N
                write(1,*) xs(i), Phi(i)
            end do

            close(1)

        end subroutine s_write_data

    ! This subroutine estimates the error via an estimate of the largest
    ! eigenvalue. Its inputs are:
    !   x1 - the solution at iteration N
    !   x2 - the solution at iteration N-1
    !   delta - a 2x1 array with the magnitudes of delta
    !   N_iter - the iteration count
    !   error - the error estimate
    subroutine s_estimate_error(x1, x2, delta, N_iter, error)

        real(kind(0d0)), dimension(0:N) :: x1, x2, diff
        real(kind(0d0)), dimension(2) :: delta
        real(kind(0d0)) :: error, lambda
        integer :: i, N_iter

        delta(2) = 0d0

        do i = 0, N
            delta(2) = delta(2) + (x1(i) - x2(i))**2d0
        end do

        delta(1) = sqrt(delta(1))

        if (N_iter > 1) then
            lambda = delta(1)/delta(2)

            error = delta(1)/(min(lambda - 1, 1d0))
        end if

        delta(1) = delta(2)

    end subroutine s_estimate_error

end module m_helpers
