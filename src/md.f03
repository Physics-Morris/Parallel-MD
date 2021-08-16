program MD
    use shared_data
    use helper
    use mpi_routines
    use mpi
    use shared_data
    use particle
    implicit none
    
    green_light = .False.

    !> checking mpi core
    call mpi_cores_check
    call mpi_comm_rank(mpi_comm_world, my_id, ierr)
    call mpi_comm_size(mpi_comm_world, numprocs, ierr)

    !> get command line argument and read input file
    !> then load particles globally
    if (my_id == master_id) then
        !> get command line argument and handle input file
        call get_cmd_arg
        !> set up particles globally, initialize particles distribution
        call print_execute_task_name('  Initialize particles distribution ... ')
        call load_particles_globally(ierr)
        call respond_to_ierr(ierr)
    end if


    !> create topology
    call print_execute_task_name('  Using Cartesian topology ... ')
    call create_cartesian_topology(ierr, numprocs)
    call respond_to_ierr(ierr)


    !> allocate particles into mpi processes
    call print_execute_task_name('  Allocating particles into mpi processes ... ')
    call allocate_particles_mpi(ierr)
    call respond_to_ierr(ierr)


    if (my_id == master_id) then
        if (green_light) then
            !> print welcome message and start the program
            call welcome
        else
            write(*, *) ' Something went wrong when initializing the program...'
            stop
        end if
    end if

    !> finishing mpi job
    call mpi_finish
end program MD
