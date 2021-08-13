program MD
    use variables
    use helper
    use mpi_routines
    use mpi
    use variables
    implicit none
    
    green_light = .False.
    !> checking mpi core
    call mpi_cores_check
    call mpi_comm_rank(mpi_comm_world, my_id, ierr)
    if (my_id .eq. master_id) then
        !> get command line argument and handle input file
        call get_cmd_arg
        !> set up particles
        call particles_initialize
        if (green_light) then
            !> print welcome message and start the program
            call welcome
        else
            write(*, *) ' Something went wrong when initializing the program...'
            stop
        end if
    end if
    call mpi_finish
end program MD
