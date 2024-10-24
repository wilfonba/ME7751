module m_vtk

    ! Dependencies ==========================================
    use m_global_parameters
    ! =======================================================

    implicit none

    private; public :: s_save_data, &
        s_open_vtk_data_file, &
        s_write_variable_to_vtk_file, &
        s_close_vtk_data_file

contains

    subroutine s_save_data(Q, n, save_count)

        type(scalar_field), dimension(1:) :: Q
        integer :: n, save_count

        call s_open_vtk_data_file(n, save_count)

        call s_write_variable_to_vtk_file(Q(1)%sf, n, 'density')
        call s_write_variable_to_vtk_file(Q(2)%sf, n, 'x-velocity')
        call s_write_variable_to_vtk_file(Q(3)%sf, n, 'y-velocity')

        call s_close_vtk_data_file()

    end subroutine s_save_data

    subroutine s_open_vtk_data_file(N, n_save)

        integer :: N, i, n_save

        character(len=100) :: file_name
        character(len=20) :: dir_name = 'data'
        character(len=100) :: line
        logical :: dir_exists

        inquire(file=trim(dir_name), exist=dir_exists)

        if (.not. dir_exists) then
            call system('mkdir '//trim(dir_name))
        end if

        file_name = trim(dir_name)//'/output_'//trim(f_int_to_str(n_save))//'.vtr'

        open (3, FILE=trim(file_name), &
              FORM='formatted', &
              STATUS='replace')

        ! header
        write(3,"(A)") "<?xml version='1.0'?>"
        write(3,"(A)") "<VTKFile type='RectilinearGrid' version='0.1' byte_order='LittleEndian'>"
        line = trim(f_int_to_str(0))//" "//trim(f_int_to_str(N))//" "// &
                trim(f_int_to_str(0))//" "//trim(f_int_to_str(N))//" "// &
                trim(f_int_to_str(0))//" "//trim(f_int_to_str(0))
        write(3,"(A)") "    <RectilinearGrid WholeExtent='"//trim(line)//"'>"
        write(3,"(A)") "        <Piece Extent='"//trim(line)//"'>"

        ! x-coordinates
        write(3,"(A)") "            <Coordinates>"
        write(3,"(A)") "                <DataArray type='Float64' format='ascii'>"
        write(3,"(A)",advance='no') "               "
        do i = 0, N
            write(3,"(A)",advance='no') trim(f_dbl_to_str(xs(i)))//" "
        end do
        write(3,*)
        write(3,"(A)") "                </DataArray>"

        ! y-coordinates
        write(3,"(A)") "                <DataArray type='Float64' format='ascii'>"
        write(3,"(A)",advance='no') "               "
        do i = 0, N
            write(3,"(A)",advance='no') trim(f_dbl_to_str(ys(i)))//" "
        end do
        write(3,*)
        write(3,"(A)") "                </DataArray>"

        ! z-coordinates
        write(3,"(A)") "                <DataArray type='Float64' format='ascii'>"
        write(3,"(A)",advance='no') "               "
        write(3,"(A)",advance='no') trim(f_dbl_to_str(0d0))//" "
        write(3,"(A)") trim(f_dbl_to_str(0d0))
        write(3,"(A)") "                </DataArray>"
        write(3,"(A)") "            </Coordinates>"

        ! point data
        write(3,"(A)") "            <PointData>"

    end subroutine s_open_vtk_data_file

    subroutine s_write_variable_to_vtk_file(Q, N, name)

        real(kind(0d0)), dimension(0:, 0:) :: Q
        integer :: i, j, N, idx
        character(len=*) :: name

        write(3,"(A)") "                <DataArray type='Float64' Name='"//trim(name)//"' format='ascii'>"
        do i = 0, N
            do j = 0, N
                write(3,"(A)",advance='no') trim(f_dbl_to_str(Q(j,i)))//" "
            end do
        end do
        write(3,*)
        write(3,"(A)") "                </DataArray>"

    end subroutine s_write_variable_to_vtk_file

    subroutine s_close_vtk_data_file()

        write(3,"(A)") "            </PointData>"
        write(3,"(A)") "        </Piece>"
        write(3,"(A)") "    </RectilinearGrid>"
        write(3,"(A)") "</VTKFile>"

        close(3)

    end subroutine s_close_vtk_data_file

    function f_int_to_str(N) result(res)

        integer, intent(in) :: N
        character(len=10) :: res

        write (res, '(I0)') N

    end function

    function f_dbl_to_str(d) result(str)

        real(8), intent(in) :: d
        character(LEN=64) :: str
        write(str, '(E16.8)') d

   end function

end module m_vtk
