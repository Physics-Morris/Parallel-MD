module mpi_routines
    use mpi
    use constants
    use helper
    use shared_data
    use particle
    use error_handle
    implicit none
    
    contains


    subroutine mpi_cores_check
        implicit none
        integer          :: numprocs, my_id, ierr
        double precision :: sq_numprocs

        call mpi_init(ierr)
        call mpi_comm_size(mpi_comm_world, numprocs, ierr)
        call mpi_comm_rank(mpi_comm_world, my_id, ierr)

        if (my_id == master_id) then
            sim_start_time = mpi_wtime()
        end if

        !> number of cpus must be a perfect square
        sq_numprocs = nint(dble(numprocs)**(1.d0/3.d0))
        if ( abs(numprocs - sq_numprocs**3) .gt. .1d0) then
            if (my_id == master_id) then
                call error_message(err_num=2)
            end if
            call mpi_finalize(ierr)
            stop
        end if
    end subroutine mpi_cores_check


    subroutine mpi_finish
        implicit none
        integer         :: my_id, ierr
        integer         :: sim_day, sim_hr, sim_min, sim_sec

        !> finishing mpi
        call mpi_comm_rank(mpi_comm_world, my_id, ierr)

        if (my_id == master_id) then
            !> unload global particle
            call print_empty_line(1)
            call print_execute_task_name('  Unloading global particles ... ', task_start_time)
            call unload_particles_globally(ierr)
            call print_task_time(task_start_time)
            call respond_to_ierr(ierr)
        end if

        if (my_id == master_id) then
            call print_execute_task_name('  Unloading local particles ... ', task_start_time)
        end if
        call unload_particles_locally(ierr)
        call print_task_time(task_start_time)
        call respond_to_ierr(ierr)

        sim_end_time = mpi_wtime()
        call mpi_finalize(ierr)
        if ((my_id == master_id).and.(ierr == 0)) then
            write(*, '(A)', advance='no') '  Finishing MPI routine ... '
            if (ierr == 0) then
                call successful
            else
                call failed
            end if
            call print_empty_line(1)
            sim_sec = int(sim_end_time-sim_start_time)
            call sec_to_day_hour_min(sim_day, sim_hr, sim_min, sim_sec)
            write(*, '(A, I0.2, A, I0.2, A, I0.2, A, I2, A5)') &
            '  Total simulation time: ', sim_hr, ':', sim_min, ':', sim_sec, '   (', sim_day, ' day)'
            call print_empty_line(1)
        end if

    end subroutine mpi_finish


    !> brocasting input file variables
    subroutine brocasting_input_file(ierr)
        implicit none
        integer, intent(out)   :: ierr
        call mpi_bcast(sim_dimension, 1, mpi_integer, master_id, mpi_comm_world, ierr)
        call mpi_bcast(x_min, 1, mpi_double_precision, master_id, mpi_comm_world, ierr)
        call mpi_bcast(x_max, 1, mpi_double_precision, master_id, mpi_comm_world, ierr)
        call mpi_bcast(y_min, 1, mpi_double_precision, master_id, mpi_comm_world, ierr)
        call mpi_bcast(y_max, 1, mpi_double_precision, master_id, mpi_comm_world, ierr)
        call mpi_bcast(z_min, 1, mpi_double_precision, master_id, mpi_comm_world, ierr)
        call mpi_bcast(z_max, 1, mpi_double_precision, master_id, mpi_comm_world, ierr)
        call mpi_bcast(total_particles, 1, mpi_integer, master_id, mpi_comm_world, ierr)
        call mpi_bcast(particle_mass, 1, mpi_double_precision, master_id, mpi_comm_world, ierr)
        call mpi_bcast(particle_charge, 1, mpi_double_precision, master_id, mpi_comm_world, ierr)
        call mpi_bcast(number_snapshots, 1, mpi_integer, master_id, mpi_comm_world, ierr)
    end subroutine brocasting_input_file
    

    !> allocate particles into mpi processes
    subroutine allocate_particles_mpi(ierr, cart_comm_3d, numprocs_x, numprocs_y, numprocs_z, local_particles)
        implicit none
        integer, intent(inout) :: ierr
        integer, intent(in)    :: cart_comm_3d
        integer, intent(in)    :: numprocs_x, numprocs_y, numprocs_z
        integer, intent(inout) :: local_particles
        integer                :: dims=3, isperiodic
        integer                :: coords(3)
        double precision       :: local_xmin, local_xmax
        double precision       :: local_ymin, local_ymax
        double precision       :: local_zmin, local_zmax

        integer                :: i, j, k
        double precision       :: pos_x, pos_y, pos_z
        integer                :: loc(3)
        integer                :: count(1:numprocs_x, 1:numprocs_y, 1:numprocs_z)
        integer                :: dest
        integer                :: status(MPI_STATUS_SIZE)
        
        !> get processor cartesian coordinate
        call mpi_cart_get(cart_comm_3d, 3, dims, isperiodic, coords, ierr)
        !> get proceccors physical min and max
        call get_local_min_max(x_min, x_max, local_xmin, local_xmax, numprocs_x, coords(1), ierr)
        call get_local_min_max(y_min, y_max, local_ymin, local_ymax, numprocs_y, coords(2), ierr)
        call get_local_min_max(z_min, z_max, local_zmin, local_zmax, numprocs_z, coords(3), ierr)

        !> let master node go through particles and send the number of particle in each procesor
        count = 0
        call mpi_comm_rank(mpi_comm_world, my_id, ierr)
        if (my_id == master_id) then
            do i = 1, total_particles
                pos_x = global_part_list(i) % pos_x
                pos_y = global_part_list(i) % pos_y
                pos_z = global_part_list(i) % pos_z
                loc = map_particle_to_global_cell(pos_x, pos_y, pos_z)
                count(loc(1), loc(2), loc(3)) = count(loc(1), loc(2), loc(3)) + 1
            end do
            !> tell each process how many particle they will recv
            dest = 0
            do i = 1, numprocs_x
                do j = 1, numprocs_y
                    do k = 1, numprocs_z
                        call mpi_send(count(i, j, k), 1, mpi_integer, dest, dest, cart_comm_3d, ierr)
                        dest = dest + 1
                    end do
                end do
            end do
        end if

        !> slaves recv number of particle
        call mpi_comm_rank(cart_comm_3d, my_id, ierr)
        call mpi_recv(local_particles, 1, mpi_integer, master_id, my_id, cart_comm_3d, status, ierr)
        call mpi_barrier(cart_comm_3d, ierr)

        !> setup local part list (allocate and send from master)
        call setup_local_part_list(local_particles, ierr)
    end subroutine allocate_particles_mpi


    !> allocate and recv local part list from master node
    subroutine setup_local_part_list(local_particles, ierr)
        implicit none
        integer, intent(in)  :: local_particles
        integer              :: my_id, ierr
        integer              :: i
        double precision     :: pos_x, pos_y, pos_z
        integer              :: loc(3), dest
        integer              :: count_recv=0

        call mpi_comm_rank(cart_comm_3d, my_id, ierr)
        !> first allocate space for storing local particles
        allocate(local_part_list(local_particles))


        !> go through all particle again and send to correspoinding processor
        call mpi_comm_rank(mpi_comm_world, my_id, ierr)
        if (my_id == master_id) then
            do i = 1, total_particles
                pos_x = global_part_list(i) % pos_x
                pos_y = global_part_list(i) % pos_y
                pos_z = global_part_list(i) % pos_z
                loc = map_particle_to_global_cell(pos_x, pos_y, pos_z)
                !> change the cell of the particle
                global_part_list(i) % global_cell_index_x = loc(1)
                global_part_list(i) % global_cell_index_y = loc(2)
                global_part_list(i) % global_cell_index_z = loc(3)
                !> this is the desitination for sending paritcle
                dest = (loc(1)-1) * numprocs_y * numprocs_z + (loc(2)-1) * numprocs_z + loc(3) - 1
                call send_particle(global_part_list(i), dest, particle_struc, ierr)
            end do
        end if

        !> receive particle
        do while (count_recv < local_particles)
            count_recv = count_recv + 1
            call recv_particle(count_recv, particle_struc, ierr)
        end do

    end subroutine setup_local_part_list


    !> create particle structure for message passing
    subroutine create_mpi_particle_struc(particle_struc, ierr)
        implicit none
        type(particle_type)                                         :: particle
        integer, intent(out)                                        :: particle_struc
        integer, intent(out)                                        :: ierr
        integer, parameter                                          :: struc_length=12
        integer, dimension(struc_length)                            :: blength
        integer(kind=mpi_address_kind), dimension(struc_length)     :: dist
        integer, dimension(struc_length+1)                          :: addr
        integer, dimension(struc_length)                            :: dtype

        integer                                                     :: i

        !> define data type contain in particle_type
        dtype(1)  = mpi_integer
        dtype(2)  = mpi_integer
        dtype(3)  = mpi_integer
        dtype(4)  = mpi_integer
        dtype(5)  = mpi_double_precision
        dtype(6)  = mpi_double_precision
        dtype(7)  = mpi_double_precision
        dtype(8)  = mpi_double_precision
        dtype(9)  = mpi_double_precision
        dtype(10) = mpi_double_precision
        dtype(11) = mpi_double_precision
        dtype(12) = mpi_double_precision

        !> block length is 1 same for all
        blength = 1

        !> get address
        call mpi_get_address(particle                      , addr(1) , ierr)
        call mpi_get_address(particle % id                 , addr(2) , ierr)
        call mpi_get_address(particle % global_cell_index_x, addr(3) , ierr)
        call mpi_get_address(particle % global_cell_index_y, addr(4) , ierr)
        call mpi_get_address(particle % global_cell_index_z, addr(5) , ierr)
        call mpi_get_address(particle % mass               , addr(6) , ierr)
        call mpi_get_address(particle % charge             , addr(7) , ierr)
        call mpi_get_address(particle % pos_x              , addr(8) , ierr)
        call mpi_get_address(particle % pos_y              , addr(9) , ierr)
        call mpi_get_address(particle % pos_z              , addr(10), ierr)
        call mpi_get_address(particle % vel_x              , addr(11), ierr)
        call mpi_get_address(particle % vel_y              , addr(12), ierr)
        call mpi_get_address(particle % vel_z              , addr(13), ierr)

        !> distance from particle
        do i = 1, struc_length
            dist(i) = addr(i+1) - addr(1)
        end do

        !> create mpi structure
        call mpi_type_create_struct(struc_length, blength, dist, dtype, particle_struc, ierr)
        !> commit type, allow mpi for optimization
        call mpi_type_commit(particle_struc, ierr)
    end subroutine create_mpi_particle_struc
    

    !> send derived particle type to other processor
    subroutine send_particle(particle, destination, particle_struc, ierr)
        implicit none
        type(particle_type), intent(in)                             :: particle
        integer, intent(in)                                         :: destination
        integer                                                     :: particle_struc
        integer, intent(out)                                        :: ierr

        call mpi_send(particle, 1, particle_struc, destination, destination, cart_comm_3d, ierr)
    end subroutine send_particle


    !> receive particle from master
    subroutine recv_particle(local_index, particle_struc, ierr)
        implicit none
        integer, intent(in)      :: local_index
        integer, intent(in)      :: particle_struc
        integer, intent(out)     :: ierr
        integer                  :: my_id
        integer                  :: status(mpi_status_size)

        call mpi_comm_rank(cart_comm_3d, my_id, ierr)
        call mpi_recv(local_part_list(local_index), 1, particle_struc, master_id, my_id, &
                      cart_comm_3d, status, ierr)
    end subroutine recv_particle


    !> create communicator 
    subroutine create_cartesian_topology(ierr, numprocs, cart_comm_3d, numprocs_x, numprocs_y, numprocs_z)
        implicit none
        integer, intent(out) :: ierr
        integer, intent(in)  :: numprocs
        integer, intent(out) :: cart_comm_3d
        integer, intent(out) :: numprocs_x, numprocs_y, numprocs_z 
        integer              :: ndims, reorder
        integer              :: dim_size(3)
        logical              :: periods(0:2)

        ndims = 3
        dim_size(1) = int(numprocs**(1.d0/3.d0))
        dim_size(2) = int(numprocs**(1.d0/3.d0))
        dim_size(3) = int(numprocs**(1.d0/3.d0))
        periods(0) = .false.
        periods(1) = .false.
        periods(2) = .false.
        reorder = 1
        !> create cartesian topology from mpi_comm_world
        call mpi_cart_create(mpi_comm_world, ndims, dim_size, &
                             periods, reorder, cart_comm_3d, ierr)
        numprocs_x = dim_size(1)
        numprocs_y = dim_size(2)
        numprocs_z = dim_size(3)
        if (my_id == master_id) then
            write(*, '(A, I1, A, I1, A, I1, A)', advance='no') &
            '(', numprocs_x, 'x', numprocs_y, 'x', numprocs_z, ')'
        end if
    end subroutine create_cartesian_topology


    !> get local minimum and maximum of simulation domain
    subroutine get_local_min_max(global_min, global_max, local_min, local_max, dim, ith, ierr)
        implicit none
        double precision, intent(in)    :: global_min, global_max
        double precision, intent(out)   :: local_min, local_max
        integer         , intent(in)    :: dim, ith
        integer         , intent(out)   :: ierr
        double precision                :: width

        if (global_min >= global_max) then
            write(*, *) global_min, global_max
            call error_message(err_num=10)
        end if
        if (ith > dim) then
            call error_message(err_num=11)
        end if

        width = (global_max - global_min) / dble(dim)
        local_min = global_min + (dble(ith) - 1.d0) * width
        local_max = global_min + dble(ith) * width
        if ((local_min >= global_min) .and. (local_max <= global_max)) then
            ierr = 0
        else
            ierr = 1
        end if
    end subroutine get_local_min_max


end module mpi_routines