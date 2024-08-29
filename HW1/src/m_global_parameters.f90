module m_global_parameters

    ! Dependencies ==========================================

    ! =======================================================

    implicit none

    ! User Inputs
    integer :: N ! number if discretization points
    real(kind(0d0)) :: gamma0 ! coeficient of zero degre gamma term
    real(kind(0d0)) :: gamma1 ! coeficient of first degree gamma term
    real(kind(0d0)) :: L ! length of domain
    real(kind(0d0)) :: Q0 ! zero order source term
    real(kind(0d0)) :: Q1 ! first order source term
    real(kind(0d0)) :: U ! Advective velocity
    real(kind(0d0)) :: tol ! Convergence tolerance for nonlinear equations

    real(kind(0d0)) :: dx ! delta x
    real(kind(0d0)), allocatable, dimension(:) :: xs ! x values

    public

    contains

        subroutine s_initialize_global_parameters()

            integer :: i ! standard loop iterator

            dx = 1d0/N

            allocate(xs(0:N))

            do i = 0,N
                xs(i) = i*dx
            end do

        end subroutine s_initialize_global_parameters

end module m_global_parameters
