module mpi_routines
    use mpi
    use constants
    use helper
    use shared_data
    use particle
    use error_handle
    use omp_lib
    implicit none

    
    contains


    subroutine mpi_setup
        implicit none
        integer          :: numprocs, my_id, ierr

        call mpi_init(ierr)
        call mpi_comm_size(mpi_comm_world, numprocs, ierr)
        call mpi_comm_rank(mpi_comm_world, my_id, ierr)

        if (my_id == master_id) then
            sim_start_time = mpi_wtime()
        end if
    end subroutine mpi_setup


    !> old routine (cores must be perfect cube)
    subroutine mpi_cube_cores_check
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
    end subroutine mpi_cube_cores_check


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
        call mpi_bcast(init_numprocs_x, 1, mpi_integer, master_id, mpi_comm_world, ierr)
        call mpi_bcast(init_numprocs_y, 1, mpi_integer, master_id, mpi_comm_world, ierr)
        call mpi_bcast(init_numprocs_z, 1, mpi_integer, master_id, mpi_comm_world, ierr)
        call mpi_bcast(load_balance, 1, mpi_logical, master_id, mpi_comm_world, ierr)
        call mpi_bcast(load_balance_num_step, 1, mpi_integer, master_id, mpi_comm_world, ierr)
        call mpi_bcast(total_particles, 1, mpi_integer, master_id, mpi_comm_world, ierr)
        call mpi_bcast(particle_mass, 1, mpi_double_precision, master_id, mpi_comm_world, ierr)
        call mpi_bcast(particle_charge, 1, mpi_double_precision, master_id, mpi_comm_world, ierr)
        call mpi_bcast(number_snapshots, 1, mpi_integer, master_id, mpi_comm_world, ierr)
        call mpi_bcast(num_auxi_per_procs, 1, mpi_integer, master_id, mpi_comm_world, ierr)
        call mpi_bcast(load_balance, 1, mpi_logical, master_id, mpi_comm_world, ierr)
        call mpi_bcast(load_balance_num_step, 1, mpi_integer, master_id, mpi_comm_world, ierr)
        call mpi_bcast(load_balance_extent, 1, mpi_double_precision, master_id, mpi_comm_world, ierr)
        call mpi_bcast(load_balance_threshold, 1, mpi_double_precision, master_id, mpi_comm_world, ierr)
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
        integer                :: status(mpi_status_size)
        
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
        integer, parameter                                          :: struc_length=13
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
        dtype(13) = mpi_integer

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
        call mpi_get_address(particle % procs_rank         , addr(14), ierr)

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

        !> finial ckeck on numprocs
        if (numprocs /= init_numprocs_x*init_numprocs_y*init_numprocs_z) then
            stop
        end if
        ndims = 3
        ! dim_size(1) = int(numprocs**(1.d0/3.d0))
        ! dim_size(2) = int(numprocs**(1.d0/3.d0))
        ! dim_size(3) = int(numprocs**(1.d0/3.d0))
        dim_size(1) = init_numprocs_x
        dim_size(2) = init_numprocs_y
        dim_size(3) = init_numprocs_z
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
            '[', numprocs_x, 'x', numprocs_y, 'x', numprocs_z, ']'
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


    !> one recv particle in local_part_list update procs rank each particle in
    subroutine update_particle_procs_rank(ierr)
        implicit none
        integer, intent(out) :: ierr
        integer              :: i, my_id

        call mpi_comm_rank(cart_comm_3d, my_id, ierr)
        do i = 1, local_particles
            local_part_list(i) % procs_rank = my_id
        end do
    end subroutine update_particle_procs_rank
    

    !> setup auxililary cell
    subroutine auxiliary_cell_setup(numprocs_x, numprocs_y, numprocs_z, ierr)
        implicit none
        integer :: numprocs_x, numprocs_y, numprocs_z
        integer :: ierr

        !> number of auxiliary cell
        auxi_num_x = numprocs_x * num_auxi_per_procs
        auxi_num_y = numprocs_y * num_auxi_per_procs
        auxi_num_z = numprocs_z * num_auxi_per_procs

        !> width of auxiliary cell
        auxi_cell_wx = (x_max - x_min) / dble(auxi_num_x)
        auxi_cell_wy = (y_max - y_min) / dble(auxi_num_y)
        auxi_cell_wz = (z_max - z_min) / dble(auxi_num_z)

        !> allocate auxiliary cell (1, 2, 3 for x, y, z; 4 for id, 5 for num part)
        allocate(auxi_cell(auxi_num_x, auxi_num_y, auxi_num_z, 5), stat=ierr)
        auxi_cell = 0
    end subroutine auxiliary_cell_setup


    !> filling up auxiliary cell with equal length at each direction by different procs
    subroutine fillup_auxi_cell(numprocs_x, numprocs_y, numprocs_z, ierr)
        implicit none
        integer          :: numprocs_x, numprocs_y, numprocs_z
        integer          :: ierr
        integer          :: i, j, k
        integer          :: procs_x, procs_y, procs_z
        double precision :: procs_wx, procs_wy, procs_wz
        double precision :: cell_cx, cell_cy, cell_cz
        integer          :: procs_id

        !> processor width in each direction (initially same width)
        procs_wx = (x_max - x_min) / dble(numprocs_x)
        procs_wy = (y_max - y_min) / dble(numprocs_y)
        procs_wz = (z_max - z_min) / dble(numprocs_z)

        do i = 1, auxi_num_x
            do j = 1, auxi_num_y
                do k = 1, auxi_num_z
                    !> center location of auxi cell
                    cell_cx = auxi_cell_wx * (dble(i) - .5d0)
                    cell_cy = auxi_cell_wy * (dble(j) - .5d0)
                    cell_cz = auxi_cell_wz * (dble(k) - .5d0)
                    procs_x = floor((cell_cx - x_min) / procs_wx)
                    procs_y = floor((cell_cy - y_min) / procs_wy)
                    procs_z = floor((cell_cz - z_min) / procs_wz)
                    procs_id = procs_x * numprocs_y * numprocs_z + procs_y * numprocs_z + procs_z
                    !> determin cell center belonging processor
                    !> (1, 2, 3) for cartesian coord of processor (start from 1)
                    !> 4 for procs_id
                    auxi_cell(i, j, k, 1) = procs_x + 1
                    auxi_cell(i, j, k, 2) = procs_y + 1
                    auxi_cell(i, j, k, 3) = procs_z + 1
                    auxi_cell(i, j, k, 4) = procs_id
                end do
            end do
        end do
        call mpi_barrier(cart_comm_3d, ierr)
    end subroutine fillup_auxi_cell


    !> allocate x, y, z slice array
    subroutine allocate_xyz_dlb_slice(ierr)
        implicit none
        integer, intent(out) :: ierr
        !> z location of the cut
        allocate(z_slice(numprocs_z-1), stat=ierr)
        !> y location of the cut for particular (z) procs
        allocate(y_slice(numprocs_z, numprocs_y-1), stat=ierr)
        !> x location of the cut for particular (z, y) procs
        allocate(x_slice(numprocs_z, numprocs_y, numprocs_x-1), stat=ierr)
        x_slice = 0
        y_slice = 0
        z_slice = 0
    end subroutine allocate_xyz_dlb_slice


    !> subroutine for dynamics load balance use rectilinear method
    !1> Collect auxiliary grid of particle number in each procs
    !2> On master node use rectilienear method for load balance
    !3> Update auxi_cell for each procs
    !4> Redistribute particle
    subroutine dynamics_load_balance(ierr, max_speedup, z_slice, y_slice, x_slice)
        implicit none
        integer, intent(out)          :: ierr
        !> number of slice needs to be made
        integer                       :: num_z_slice, num_y_slice, num_x_slice
        !> z location of the cut
        integer, intent(out)          :: z_slice(numprocs_z-1)
        !> y location of the cut for particular (z) procs
        integer, intent(out)          :: y_slice(numprocs_z, numprocs_y-1)
        !> x location of the cut for particular (z, y) procs
        integer, intent(out)          :: x_slice(numprocs_z, numprocs_y, numprocs_x-1)
        integer                       :: total_auxi
        integer                       :: before_max, after_max
        double precision, intent(out) :: max_speedup

        total_auxi = auxi_num_x*auxi_num_y*auxi_num_z

        allocate(auxi_cell_new(auxi_num_x, auxi_num_y, auxi_num_z, 5), stat=ierr)

        num_z_slice = numprocs_z-1
        num_y_slice = numprocs_y-1
        num_x_slice = numprocs_x-1

        !> get current processor id
        call mpi_comm_rank(cart_comm_3d, my_id, ierr)

        !1> Collect auxiliary grid of particle number in each procs
        call get_auxi_cell_part_num(ierr)

        if (my_id == master_id) then
            !2> On master node use rectilienear method for load balance
            !2> get the cut in each direction
            call get_dlb_slice_xyz(z_slice, y_slice, x_slice)

            !3> create new auxi cell accroding to the slice we found
            call update_auxi_cell(z_slice, y_slice, x_slice, ierr)
        end if

        !4> Redistribute particle
        auxi_cell_new(:, :, :, 5) = auxi_cell(:, :, :, 5)
        call mpi_bcast(auxi_cell_new(:, :, :, :), total_auxi*5, &
                       mpi_integer, master_id, cart_comm_3d, ierr)
        call dlb_redistribute_particles(ierr, before_max, after_max)

        !> calculate maximum speedup
        max_speedup = dble(before_max) / dble(after_max)

        deallocate(auxi_cell_new, stat=ierr)
    end subroutine dynamics_load_balance


    !> collect auxi cell particle number
    subroutine get_auxi_cell_part_num(ierr)
        integer               :: i
        integer, intent(out)  :: ierr
        integer, allocatable  :: local_num(:, :, :)
        double precision      :: x, y, z
        integer               :: cell_x, cell_y, cell_z
        integer               :: total_auxi

        call mpi_comm_rank(cart_comm_3d, my_id, ierr)
        allocate(local_num(auxi_num_x, auxi_num_y, auxi_num_z))
        local_num = 0

        !> first loop over local particle and fill up local_num
        do i = 1, local_particles
            x = local_part_list(i) % pos_x
            y = local_part_list(i) % pos_y
            z = local_part_list(i) % pos_z
            !> find out which auxi cell part in
            cell_x = floor((x - x_min) / auxi_cell_wx) + 1
            cell_y = floor((y - y_min) / auxi_cell_wy) + 1
            cell_z = floor((z - z_min) / auxi_cell_wz) + 1
            local_num(cell_x, cell_y, cell_z) = local_num(cell_x, cell_y, cell_z) + 1
        end do

        !> after caluclate locally master collect all the data from procs in to auxi_cell
        auxi_cell(:, :, :, 5) = 0
        total_auxi = auxi_num_x * auxi_num_y * auxi_num_z
        call mpi_allreduce(local_num, auxi_cell(:, :, :, 5), total_auxi, mpi_integer, &
                           mpi_sum, cart_comm_3d, ierr)
        deallocate(local_num)
    end subroutine get_auxi_cell_part_num
    

    !> get slice for DLB
    subroutine get_dlb_slice_xyz(z_slice, y_slice, x_slice)
        implicit none
        integer               :: i, j, k
        !> number of slice needs to be made
        integer               :: num_z_slice, num_y_slice, num_x_slice
        !> z location of the cut
        integer, intent(out)  :: z_slice(numprocs_z-1)
        !> y location of the cut for particular (z) procs
        integer, intent(out)  :: y_slice(numprocs_z, numprocs_y-1)
        !> x location of the cut for particular (z, y) procs
        integer, intent(out)  :: x_slice(numprocs_z, numprocs_y, numprocs_x-1)
        !> target number of the cut
        integer               :: target_part
        !> keep count of current particle
        integer               :: current_part
        !> keep count of the number of slice created
        integer               :: current_slice
        !> the range in z direction of the start and the end of the domain
        integer               :: z_dir_start, z_dir_end
        !> the range in y direction of the start and the end of the domain
        integer               :: y_dir_start, y_dir_end

        num_z_slice = numprocs_z-1
        num_y_slice = numprocs_y-1
        num_x_slice = numprocs_x-1

        z_slice = 0
        y_slice = 0
        x_slice = 0

        !> start from z-direction (create numprocs_z-1 slice)
        if (num_z_slice /= 0) then
            target_part = sum(auxi_cell(:, :, :, 5)) / numprocs_z
            current_part = 0
            current_slice = 1
            !> go through z direction auxiliary grid to find slice
            do i = 1, auxi_num_z 
                !> add a slice of x-y plane particle(s)
                current_part = current_part + sum(auxi_cell(:, :, i, 5))
                !> until it reach desire value
                if (current_part >= target_part) then
                    !> rememeber the location of that slice
                    z_slice(current_slice) = i
                    current_slice = current_slice + 1
                    !> if there were enough slice exit, else next slice
                    if (current_slice > num_z_slice) then
                        exit
                    end if
                    !> reset particle number and continue search for next slice
                    current_part = 0
                end if
            end do
        end if

        !> then in every z-slice x-y retangle, create (num_y_slice) slice in y direction
        if (num_y_slice /= 0) then
            do i = 1, numprocs_z
                !> start and end of the slice in z [start, end)
                if (i == 1) then
                    z_dir_start = 1
                else
                    z_dir_start = z_slice(i-1)
                end if
                if (i == numprocs_z) then
                    z_dir_end = auxi_num_z
                else
                    z_dir_end = z_slice(i) - 1
                end if

                target_part = sum(auxi_cell(:, :, z_dir_start:z_dir_end, 5)) / numprocs_y
                current_part = 0
                current_slice = 1
                !> go through y direction auxiliary grid to find slice
                do j = 1, auxi_num_y 
                    !> add a slice of x-y plane particle(s)
                    current_part = current_part + sum(auxi_cell(:, j, z_dir_start:z_dir_end, 5))
                    !> until it reach desire value
                    if (current_part >= target_part) then
                        !> rememeber the location of that slice
                        y_slice(i, current_slice) = j
                        current_slice = current_slice + 1
                        !> if there were enough slice exit, else next slice
                        if (current_slice > num_y_slice) then
                            exit
                        end if
                        !> reset particle number and continue search for next slice
                        current_part = 0
                    end if
                end do
            end do
        end if

        !> finally x-direction (in every z,y slice create numprocs_x-1 slice)
        if (num_x_slice /= 0) then
            do i = 1, numprocs_z
                if (i == 1) then
                    z_dir_start = 1
                else
                    z_dir_start = z_slice(i-1)
                end if
                if (i == numprocs_z) then
                    z_dir_end = auxi_num_z
                else
                    z_dir_end = z_slice(i) - 1
                end if
                do j = 1, numprocs_y
                    !> start and end of the slice in z [start, end)
                    !> start and end of the slice in y [start, end)
                    if (j == 1) then
                        y_dir_start = 1
                    else
                        y_dir_start = y_slice(i, j-1)
                    end if
                    if (j == numprocs_y) then
                        y_dir_end = auxi_num_y
                    else
                        y_dir_end = y_slice(i, j) - 1
                    end if

                    target_part = sum(auxi_cell(:, y_dir_start:y_dir_end, z_dir_start:z_dir_end, 5)) &
                                    / numprocs_x
                    current_part = 0
                    current_slice = 1
                    !> go through y direction auxiliary grid to find slice
                    do k = 1, auxi_num_x
                        !> add a slice of x-y plane particle(s)
                        current_part = current_part + &
                                        sum(auxi_cell(k, y_dir_start:y_dir_end, z_dir_start:z_dir_end, 5))
                        !> until it reach desire value
                        if (current_part >= target_part) then
                            !> rememeber the location of that slice
                            x_slice(i, j, current_slice) = k
                            current_slice = current_slice + 1
                            !> if there were enough slice exit, else next slice
                            if (current_slice > num_x_slice) then
                                exit
                            end if
                            !> reset particle number and continue search for next slice
                            current_part = 0
                        end if
                    end do
                end do
            end do
        end if
    
    end subroutine get_dlb_slice_xyz


    !> create new auxi cell accroading to the dle slice in x, y, z
    subroutine update_auxi_cell(z_slice, y_slice, x_slice, ierr)
        implicit none
        integer, intent(out)        :: ierr
        !> z location of the cut
        integer, intent(in)         :: z_slice(numprocs_z-1)
        !> y location of the cut for particular (z) procs
        integer, intent(in)         :: y_slice(numprocs_z, numprocs_y-1)
        !> x location of the cut for particular (z, y) procs
        integer, intent(in)         :: x_slice(numprocs_z, numprocs_y, numprocs_x-1)
        integer                     :: z_start, z_end, y_start, y_end, x_start, x_end
        integer                     :: procs_z, procs_y, procs_x
        integer                     :: procs_id

        
        do procs_z = 1, numprocs_z
            !> get start and end in z direction
            if (procs_z == 1) then
                z_start = 1
            else
                z_start = z_slice(procs_z-1)
            end if
            if (procs_z == numprocs_z) then
                z_end = auxi_num_z
            else
                z_end = z_slice(procs_z) - 1
            end if
            do procs_y = 1, numprocs_y
                !> get start and end in y direction
                if (procs_y == 1) then
                    y_start = 1
                else
                    y_start = y_slice(procs_z, procs_y-1)
                end if
                if (procs_y == numprocs_y) then
                    y_end = auxi_num_y
                else
                    y_end = y_slice(procs_z, procs_y) - 1
                end if
                do procs_x = 1, numprocs_x
                    !> get start and end in x direction
                    if (procs_x == 1) then
                        x_start = 1
                    else
                        x_start = x_slice(procs_z, procs_y, procs_x-1)
                    end if
                    if (procs_x == numprocs_x) then
                        x_end = auxi_num_x
                    else
                        x_end = x_slice(procs_z, procs_y, procs_x) - 1
                    end if

                    procs_id = (procs_x-1) * numprocs_y * numprocs_z + (procs_y-1) * numprocs_z + (procs_z-1)
                    auxi_cell_new(x_start:x_end, y_start:y_end, z_start:z_end, 1) = procs_x
                    auxi_cell_new(x_start:x_end, y_start:y_end, z_start:z_end, 2) = procs_y
                    auxi_cell_new(x_start:x_end, y_start:y_end, z_start:z_end, 3) = procs_z
                    auxi_cell_new(x_start:x_end, y_start:y_end, z_start:z_end, 4) = procs_id
                end do
            end do
        end do
        ierr = 0
    end subroutine update_auxi_cell


    !> redistribute particles accroding to auxi_cell_new
    subroutine dlb_redistribute_particles(ierr, before_max, after_max)
        implicit none
        integer, intent(out)  :: ierr
        integer               :: local_particles_new
        integer               :: i, j, k
        integer               :: part_auxi(3), part_not_moving
        integer               :: send_part, recv_part
        integer               :: num_send, num_recv
        integer               :: dest
        integer, allocatable  :: status(:, :)
        integer, allocatable  :: requests(:)
        !> before and after particle maximum number
        integer, intent(out)  :: before_max, after_max

        call mpi_comm_rank(cart_comm_3d, my_id, ierr)

        !> count new local particle number
        local_particles_new = 0
        do i = 1, auxi_num_x
            do j = 1, auxi_num_y
                do k = 1, auxi_num_z
                    if (auxi_cell_new(i, j, k, 4) == my_id) then
                        local_particles_new = local_particles_new + &
                                              auxi_cell_new(i, j, k, 5)
                    end if
                end do
            end do
        end do

        !> create enough space for all particle
        allocate(local_part_list_new(local_particles_new))

        !> put the particle that doesn't need to redsitribute into the list
        part_not_moving = 0
        do i = 1, local_particles
            part_auxi = which_auxi_cell(local_part_list(i))
            dest = auxi_cell_new(part_auxi(1), part_auxi(2), part_auxi(3), 4) 
            if (dest == my_id) then
                part_not_moving = part_not_moving + 1
                local_part_list_new(part_not_moving) = local_part_list(i)
            end if
        end do

        !> start to move the particles to desitination
        send_part = local_particles - part_not_moving
        recv_part = local_particles_new - part_not_moving
        num_send = 0

        allocate(status(mpi_status_size, recv_part))
        allocate(requests(recv_part))


        do num_recv = 1, recv_part
            call mpi_irecv(local_part_list_new(part_not_moving+num_recv), 1, particle_struc, &
                           mpi_any_source, mpi_any_tag, cart_comm_3d, requests(num_recv), ierr)
        end do
        do i = 1, local_particles
            part_auxi = which_auxi_cell(local_part_list(i))
            dest = auxi_cell_new(part_auxi(1), part_auxi(2), part_auxi(3), 4) 
            if ( dest /= my_id) then
                call mpi_send(local_part_list(i), 1, particle_struc, dest, &
                              my_id, cart_comm_3d, ierr)
            end if
        end do
        call mpi_waitall(recv_part, requests, status, ierr)

        !> find before max and after max
        call mpi_allreduce(local_particles, before_max, 1, mpi_integer, mpi_max, cart_comm_3d, ierr)
        call mpi_allreduce(local_particles_new, after_max, 1, mpi_integer, mpi_max, cart_comm_3d, ierr)

        !> update new local particles number
        local_particles = local_particles_new

        !> transfer new local particle list
        deallocate(local_part_list)
        allocate(local_part_list(local_particles))
        local_part_list = local_part_list_new
        
        !> and free temporary list
        deallocate(local_part_list_new)

    end subroutine dlb_redistribute_particles


    !> check load balance threshold
    subroutine check_dlb_threshold(threshold)
        implicit none
        double precision, intent(out) :: threshold
        integer                       :: max_part, min_part
        call mpi_allreduce(local_particles, max_part, 1, mpi_integer, mpi_max, cart_comm_3d, ierr)
        call mpi_allreduce(local_particles, min_part, 1, mpi_integer, mpi_min, cart_comm_3d, ierr)
        threshold = dble(min_part) / dble(max_part)
    end subroutine check_dlb_threshold


end module mpi_routines