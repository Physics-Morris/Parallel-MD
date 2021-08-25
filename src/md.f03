program MD
    use shared_data
    use helper
    use mpi_routines
    use mpi
    use shared_data
    use particle
    use error_handle
    use diagnostics
    implicit none
    
    green_light = .False.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> checking mpi core
    call mpi_cores_check
    call mpi_comm_rank(mpi_comm_world, my_id, ierr)
    call mpi_comm_size(mpi_comm_world, numprocs, ierr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> get command line argument and read input file
    !> then load particles globally
    if (my_id == master_id) then
        !> get command line argument and handle input file
        call get_cmd_arg
        !> set up particles globally, initialize particles distribution
        call print_execute_task_name('  Initialize particles distribution ... ', task_start_time)
        call load_particles_globally(ierr)
        call print_task_time(task_start_time)
        call respond_to_ierr(ierr)
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> brocasting input file shared variables
    call print_execute_task_name('  Brocasting input file variables ... ', task_start_time)
    call brocasting_input_file(ierr)
    call print_task_time(task_start_time)
    call respond_to_ierr(ierr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> create topology
    call print_execute_task_name('  Using Cartesian topology ... ', task_start_time)
    call create_cartesian_topology(ierr, numprocs, cart_comm_3d, numprocs_x, numprocs_y, numprocs_z)
    call print_task_time(task_start_time)
    call respond_to_ierr(ierr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> create particle struc for mpi send recv
    call print_execute_task_name('  Submitting particle structure ... ', task_start_time)
    call create_mpi_particle_struc(particle_struc, ierr)
    call print_task_time(task_start_time)
    call respond_to_ierr(ierr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> allocate particles into mpi processes
    call print_execute_task_name('  Allocating particles into mpi processes ... ', task_start_time)
    call allocate_particles_mpi(ierr, cart_comm_3d, numprocs_x, numprocs_y, numprocs_z, local_particles)
    call print_task_time(task_start_time)
    call respond_to_ierr(ierr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> initial output
    call print_execute_task_name('  Write initial setup to H5 file ... ', task_start_time)
    call output(step, ierr)
    call print_task_time(task_start_time)
    call respond_to_ierr(ierr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> print welcome message and wait for all processor to join
    if (my_id == master_id) then
        if (green_light) then
            !> print welcome message and start the program
            call welcome
        else
            call error_message(err_num=1)
        end if
    end if
    call mpi_barrier(cart_comm_3d, ierr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> finishing
    call mpi_finish
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program MD