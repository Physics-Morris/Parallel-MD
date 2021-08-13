module mpi_routines
    use mpi
    use constants
    use helper
    use variables
    implicit none
    
    contains


    subroutine mpi_cores_check
        implicit none
        integer          :: numprocs, my_id, ierr
        double precision :: sq_numprocs

        call mpi_init(ierr)
        call mpi_comm_size(mpi_comm_world, numprocs, ierr)
        call mpi_comm_rank(mpi_comm_world, my_id, ierr)

        !> number of cpus must be a perfect square
        sq_numprocs = int(dsqrt(dble(numprocs)))
        if ( abs(numprocs - sq_numprocs**2) .gt. 0.1) then
            if (my_id == master_id) then
                CALL set_term_color(3)
                call print_empty_line(1)
                write(*, *)' The number of cores must be a perfect square', &
                           ' (e.g. 1, 4, 9, 16, etc...)'
                CALL set_term_color(12)
                call print_empty_line(1)
                write(*, *) ' Finishing the program ... '
                call print_empty_line(1)
            end if
            call mpi_finalize(ierr)
            stop
        end if
    end subroutine mpi_cores_check


    subroutine mpi_finish
        implicit none
        integer         :: my_id, ierr

        call mpi_comm_rank(mpi_comm_world, my_id, ierr)
        if (my_id == master_id) then
            call set_term_color(3)
            call print_empty_line(1)
            write(*, '(A)', advance='no') '  Finishing MPI routine ... '
            call successful
            call set_term_color(12)
            call print_empty_line(1)
        end if
        call mpi_finalize(ierr)
    end subroutine mpi_finish


    subroutine particles_initialize
        implicit none
        integer         :: numprocs, ierr

        call mpi_comm_size(mpi_comm_world, numprocs, ierr)
        write(*, '(A)', advance='no') '  Allocating particles into processes ... '
        call successful
        ! call failed
    end subroutine particles_initialize


end module mpi_routines