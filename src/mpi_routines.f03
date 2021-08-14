module mpi_routines
    use mpi
    use constants
    use helper
    use shared_data
    use particle
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

        !> finishing mpi
        call mpi_comm_rank(mpi_comm_world, my_id, ierr)

        if (my_id == master_id) then
            !> unload global particle
            call print_empty_line(1)
            write(*, '(A)', advance='no') '  Unloading gloabl particles ... '
            call unload_particles_globally(ierr)
            call respond_to_ierr(ierr)
        end if

        call mpi_finalize(ierr)
        if ((my_id == master_id).and.(ierr == 0)) then
            call print_empty_line(1)
            write(*, '(A)', advance='no') '  Finishing MPI routine ... '
            call respond_to_ierr(ierr)
            call print_empty_line(1)
        end if

    end subroutine mpi_finish


    subroutine particles_initialize
        implicit none
        integer         :: numprocs, ierr

        call mpi_comm_size(mpi_comm_world, numprocs, ierr)

        !> initialize particles distribution
        call load_particles_globally(ierr)
        write(*, '(A)', advance='no') '  Initialize particles distribution ... '
        call respond_to_ierr(ierr)

        !> allocate particles into mpi processes
        write(*, '(A)', advance='no') '  Allocating particles into mpi processes ... '
        call allocate_particles_mpi(ierr)
        call respond_to_ierr(ierr)
    end subroutine particles_initialize
    

    !> allocate particles into mpi processes
    subroutine allocate_particles_mpi(ierr)
        implicit none
        integer, intent(inout) :: ierr
        ierr = 0
    end subroutine allocate_particles_mpi


end module mpi_routines