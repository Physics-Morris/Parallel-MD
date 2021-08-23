module error_handle
    use mpi
    use constants
    use shared_data
    implicit none

    contains

        subroutine error_message(err_num)
            implicit none
            integer, intent(in) :: err_num
            integer             :: my_id, ierr

            call mpi_comm_rank(mpi_comm_world, my_id, ierr)

            select case (err_num)
                case(1)
                    write(*, *) ' error 1: Something went wrong when initializing the program...'
                    stop
                case(2)
                    write(*, *)
                    write(*, *)' error 2: The number of cores must be a perfect cube', &
                               ' (e.g. 1, 8, 27, 64, etc...)'
                    write(*, *)
                    write(*, *) ' Finishing the program ... '
                    write(*, *)
                case(10)
                    write(*, *)
                    write(*, *) ' error 10: When calling subroutine get_local_min_max at proccesor id: ', my_id
                    stop
                case(11)
                    write(*, *)
                    write(*, *) ' error 11: When calling subroutine get_local_min_max at proccesor id: ', my_id
                    stop
            end select

        end subroutine error_message

end module error_handle