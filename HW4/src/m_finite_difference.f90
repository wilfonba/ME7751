module m_finite_difference

    ! Dependencies ==========================================
    use m_global_parameters

    use m_helpers
    ! =======================================================

    private; public :: s_compute_rhs, &
        s_compute_lhs

contains

    subroutine s_compute_rhs(Q, RHS, N, h)

        type(scalar_field), dimension(1:) :: Q, RHS
        integer :: i, j, N
        real(kind(0d0)), dimension(-1:1) :: neighbors
        real(kind(0d0)), dimension(-1:1) :: pressures
        real(kind(0d00)) :: h

        do j = 1, N-1
            do i = 1, N-1
                ! Density equation ============================================
                RHS(1)%sf(i,j) = (dt/(beta*2d0*h)) * ( &
                    Q(2)%sf(i+1,j) - Q(2)%sf(i-1,j) + &
                    Q(2)%sf(i,j+1) - Q(2)%sf(i,j-1))

                ! x-velocity equation x-derivatives ===========================
                pressures(-1) = Q(1)%sf(i-1,j)/beta
                pressures(1)  = Q(1)%sf(i+1,j)/beta

                neighbors(-1) = Q(2)%sf(i-1,j)**2d0 + pressures(-1)
                neighbors(1) = Q(2)%sf(i+1,j)**2d0 + pressures(1)

                RHS(2)%sf(i,j) = (dt/(2d0*h)) * ( &
                    neighbors(1) - neighbors(-1))

                ! x-velocity equation y-derivatives
                neighbors(-1) = Q(2)%sf(i,j-1)*Q(3)%sf(i,j-1)
                neighbors(1) = Q(2)%sf(i,j+1)*Q(3)%sf(i, j+1)

                RHS(2)%sf(i,j) = RHS(2)%sf(i,j) + (dt/(2d0*h)) * ( &
                    neighbors(1) - neighbors(-1))

                ! y-velocity equation x-derivatives ===========================
                neighbors(-1) = Q(2)%sf(i,j-1)*Q(3)%sf(i,j-1)
                neighbors(1) = Q(2)%sf(i,j+1)*Q(3)%sf(i,j+1)

                RHS(3)%sf(i,j) = (dt/(2d0*h)) * ( &
                    neighbors(1) - neighbors(-1))

                ! y-velocity equation y-derivatives ===========================
                pressures(-1) = Q(1)%sf(i,j-1)/beta
                pressures(1)  = Q(1)%sf(i,j+1)/beta

                neighbors(-1) = Q(3)%sf(i,j-1)**2d0 + pressures(-1)
                neighbors(1) = Q(3)%sf(i,j+1)**2d0 + pressures(1)

                RHS(3)%sf(i,j) = RHS(3)%sf(i,j) + (dt/(2d0*h)) * ( &
                    neighbors(1) - neighbors(-1))

            end do
        end do

    end subroutine s_compute_rhs

    subroutine s_compute_lhs(Q, LHS, dir, N, h)

        type(scalar_field), dimension(1:) :: Q, LHS
        integer :: dir, i, j, N
        real(kind(0d0)) :: h

        real(kind(0d0)), dimension(-1:1) :: neighbors, pressures

        if (dir == 1) then
            do j = 1, N-1
                do i = 1, N-1

                    LHS(1)%sf(i,j) = 1d0
                    LHs(2)%sf(i,j) = 1d0
                    LHS(3)%sf(i,j) = 1d0

                    ! A-term density equation =================================
                    neighbors(-1) = Q(2)%sf(i-1,j)/beta
                    neighbors(1)  = Q(2)%sf(i+1,j)/beta

                    LHS(1)%sf(i,j) = LHS(1)%sf(i,j) + (dt/(4d0*h*beta)) * &
                        (neighbors(1) - neighbors(-1))

                    ! A-term x-velocity =======================================
                    pressures(-1) = Q(1)%sf(i-1,j)/beta
                    pressures(1)  = Q(1)%sf(i+1,j)/beta

                    neighbors(-1) = pressures(-1) + 2d0*Q(2)%sf(i-1,j)**2d0
                    neighbors(1)  = pressures(1) + 2d0*Q(2)%sf(i+1,j)**2d0

                    LHS(2)%sf(i,j) = LHS(2)%sf(i,j) + (dt/(4d0*h)) * &
                        (neighbors(1) - neighbors(-1))

                    ! A-term y-velocity =======================================
                    neighbors(-1) = Q(2)%sf(i-1,j)*Q(3)%sf(i-1,j)
                    neighbors(1)  = Q(2)%sf(i+1,j)*Q(3)%sf(i+1,j)

                    LHS(3)%sf(i,j) = LHS(2)%sf(i,j) + (dt/(4d0*h)) * &
                        (neighbors(1) - neighbors(-1))

                    ! D-term x-velocity =======================================
                    LHS(2)%sf(i,j) = LHS(2)%sf(i,j) + (dt/(4*Re*h**2d0)) * ( &
                        LHS(2)%sf(i+1,j) - 2d0*LHS(2)%sf(i,j) + LHS(2)%sf(i-1,j))

                    ! D-term y-velocity =======================================
                    LHS(3)%sf(i,j) = LHS(3)%sf(i,j) + (dt/(4*Re*h**2d0)) * ( &
                        LHS(3)%sf(i+1,j) - 2d0*LHS(3)%sf(i,j) + LHS(3)%sf(i-1,j))

                end do
            end do

        elseif (dir == 2) then

            do j = 1, N-1
                do i = 1, N-1

                    LHS(1)%sf(i,j) = 1d0
                    LHs(2)%sf(i,j) = 1d0
                    LHS(3)%sf(i,j) = 1d0

                    ! A-term density equation =================================
                    neighbors(-1) = Q(3)%sf(i,j-1)/beta
                    neighbors(1)  = Q(3)%sf(i,j+1)/beta

                    LHS(1)%sf(i,j) = LHS(1)%sf(i,j) + (dt/(4d0*h*beta)) * &
                        (neighbors(1) - neighbors(-1))

                    ! A-term x-velocity =======================================
                    neighbors(-1) = 2d0*Q(2)%sf(i,j-1)*Q(3)%sf(i,j-1)
                    neighbors(1)  = 2d0*Q(2)%sf(i,j+1)*Q(3)%sf(i,j+1)

                    LHS(2)%sf(i,j) = LHS(2)%sf(i,j) + (dt/(4d0*h)) * &
                        (neighbors(1) - neighbors(-1))

                    ! A-term y-velocity =======================================
                    pressures(-1) = Q(1)%sf(i,j-1)/beta
                    pressures(1)  = Q(1)%sf(i,j+1)/beta

                    neighbors(-1) = pressures(-1) + 2d0*Q(3)%sf(i,j-1)**2d0
                    neighbors(1)  = pressures(1) + 2d0*Q(3)%sf(i,j+1)**2d0

                    LHS(3)%sf(i,j) = LHS(2)%sf(i,j) + (dt/(4d0*h)) * &
                        (neighbors(1) - neighbors(-1))

                    ! D-term x-velocity =======================================
                    LHS(2)%sf(i,j) = LHS(2)%sf(i,j) + (dt/(4*Re*h**2d0)) * ( &
                        LHS(2)%sf(i,j+1) - 2d0*LHS(2)%sf(i,j) + LHS(2)%sf(i,j-1))

                    ! D-term y-velocity =======================================
                    LHS(3)%sf(i,j) = LHS(3)%sf(i,j) + (dt/(4*Re*h**2d0)) * ( &
                        LHS(3)%sf(i,j+1) - 2d0*LHS(3)%sf(i,j) + LHS(3)%sf(i,j-1))

                end do
            end do

        end if

    end subroutine s_compute_lhs

end module m_finite_difference
