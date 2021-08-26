module diagnostics
    use hdf5
    use mpi
    use shared_data
    use constants
    use math
    use error_handle
    use helper
    use particle

    implicit none

    contains


    !> write particle information to h5 file
    subroutine output(num, ierr)
        implicit none
        integer, intent(in)                                    :: num
        integer, intent(out)                                   :: ierr
        character(len=10)                                      :: file_number
        character(len=8)                                       :: fmt='(I4.4)'
        character(len=8)                                       :: default_dir='../data/'
        character(len=16)                                      :: filename
        !> dataset name
        character(len=11), parameter                           :: dsetname_id = "particle_id" 
        character(len=13), parameter                           :: dsetname_mass = "particle_mass" 
        character(len=15), parameter                           :: dsetname_charge = "particle_charge" 
        character(len=17), parameter                           :: dsetname_vel = "particle_velocity" 
        character(len=17), parameter                           :: dsetname_pos = "particle_position" 
        character(len=17), parameter                           :: dsetname_index = "global_cell_index" 
        !> file identifier 
        integer(hid_t)                                         :: file_id       
        !> dataset identifier 
        integer(hid_t)                                         :: dsetid_id
        integer(hid_t)                                         :: dsetid_mass
        integer(hid_t)                                         :: dsetid_charge
        integer(hid_t)                                         :: dsetid_vel
        integer(hid_t)                                         :: dsetid_pos
        integer(hid_t)                                         :: dsetid_index
        !> dataspace identifier in file 
        integer(hid_t)                                         :: filespace_id     
        integer(hid_t)                                         :: filespace_mass
        integer(hid_t)                                         :: filespace_charge
        integer(hid_t)                                         :: filespace_vel     
        integer(hid_t)                                         :: filespace_pos     
        integer(hid_t)                                         :: filespace_index     
        !> dataspace identifier in memory
        integer(hid_t)                                         :: memspace_id
        integer(hid_t)                                         :: memspace_mass
        integer(hid_t)                                         :: memspace_charge
        integer(hid_t)                                         :: memspace_vel
        integer(hid_t)                                         :: memspace_pos
        integer(hid_t)                                         :: memspace_index
        !> property list identifier 
        integer(hid_t)                                         :: plist_id      
        !> dataset dimensions.
        integer(hsize_t), dimension(2)                         :: dim_id
        integer(hsize_t), dimension(2)                         :: dim_mass
        integer(hsize_t), dimension(2)                         :: dim_charge
        integer(hsize_t), dimension(2)                         :: dim_vel
        integer(hsize_t), dimension(2)                         :: dim_pos
        integer(hsize_t), dimension(2)                         :: dim_index

        integer(hsize_t), dimension(2)                         :: count_id
        integer(hsize_t), dimension(2)                         :: count_mass
        integer(hsize_t), dimension(2)                         :: count_charge
        integer(hsize_t), dimension(2)                         :: count_vel
        integer(hsize_t), dimension(2)                         :: count_pos
        integer(hsize_t), dimension(2)                         :: count_index

        integer(hssize_t), dimension(2)                        :: offset_id
        integer(hssize_t), dimension(2)                        :: offset_mass
        integer(hssize_t), dimension(2)                        :: offset_charge
        integer(hssize_t), dimension(2)                        :: offset_vel
        integer(hssize_t), dimension(2)                        :: offset_pos
        integer(hssize_t), dimension(2)                        :: offset_index

        !> data to write
        double precision, allocatable                          :: data_id(:, :)
        double precision, allocatable                          :: data_mass(:, :)
        double precision, allocatable                          :: data_charge(:, :)
        double precision, allocatable                          :: data_vel(:, :)
        double precision, allocatable                          :: data_pos(:, :)
        double precision, allocatable                          :: data_index(:, :)

        !> dataset rank 
        integer                                                :: rank_id = 2
        integer                                                :: rank_mass = 2
        integer                                                :: rank_charge = 2
        integer                                                :: rank_vel = 2
        integer                                                :: rank_pos = 2
        integer                                                :: rank_index = 2

        !> error flags
        integer :: error

        !> mpi related
        integer                                                :: comm, info, dims, isperiodic, coords(3)

        !> local particle and global particle
        integer, dimension(numprocs_x, numprocs_y, numprocs_z) :: locals_part_num
        integer                                                :: global_part_num
        integer                                                :: count_part
        integer, dimension(numprocs)                           :: counts_part
        integer                                                :: part_offset

        integer                                                :: i

        !> first call mpi related subroutine
        comm = mpi_comm_world
        info = mpi_info_null
        call mpi_comm_size(comm, numprocs, ierr)
        call mpi_comm_rank(comm, my_id, ierr)
        call mpi_cart_get(cart_comm_3d, 3, dims, isperiodic, coords, ierr)

        !> output file name
        write(file_number, fmt) num
        filename = default_dir // trim(file_number) // '.h5'

        !> get all the local particle number and global particle number
        locals_part_num = count_local_particle(ierr)
        global_part_num = count_global_particle(ierr)

        !> gather and bcast particle in each processor to all processor (allgather)
        count_part = locals_part_num(coords(1)+1, coords(2)+1, coords(3)+1)
        call mpi_allgather(count_part, 1, mpi_integer, counts_part, 1, mpi_integer, mpi_comm_world, ierr)

        !> calculate particle offset with different process
        part_offset = calculate_offset(counts_part, my_id)

        !> dataset dimension
        dim_id(1)     = 1
        dim_id(2)     = global_part_num
        dim_mass(1)   = 1
        dim_mass(2)   = global_part_num
        dim_charge(1) = 1
        dim_charge(2) = global_part_num
        dim_vel(1)    = 3
        dim_vel(2)    = global_part_num
        dim_pos(1)    = 3
        dim_pos(2)    = global_part_num
        dim_index(1)  = 3
        dim_index(2)  = global_part_num

        !> initialize fortran predefined datatypes
        call h5open_f(error) 

        !> setup file access property list with parallel i/o access.
        call h5pcreate_f(h5p_file_access_f, plist_id, error)
        call h5pset_fapl_mpio_f(plist_id, comm, mpi_info_null, error)

        !> create the file collectively.
        call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error, access_prp=plist_id)
        call h5pclose_f(plist_id, error)

        !> create the data space for the dataset. 
        call h5screate_simple_f(rank_id,     dim_id,     filespace_id,     error)
        call h5screate_simple_f(rank_mass,   dim_mass,   filespace_mass,   error)
        call h5screate_simple_f(rank_charge, dim_charge, filespace_charge, error)
        call h5screate_simple_f(rank_vel,    dim_vel,    filespace_vel,    error)
        call h5screate_simple_f(rank_pos,    dim_pos,    filespace_pos,    error)
        call h5screate_simple_f(rank_index,  dim_index,  filespace_index,  error)

        !> create the dataset with default properties.
        call h5dcreate_f(file_id, dsetname_id,     h5t_ieee_f64le, filespace_id,     dsetid_id,     error)
        call h5dcreate_f(file_id, dsetname_mass,   h5t_ieee_f64le, filespace_mass,   dsetid_mass,   error)
        call h5dcreate_f(file_id, dsetname_charge, h5t_ieee_f64le, filespace_charge, dsetid_charge, error)
        call h5dcreate_f(file_id, dsetname_vel,    h5t_ieee_f64le, filespace_vel,    dsetid_vel,    error)
        call h5dcreate_f(file_id, dsetname_pos,    h5t_ieee_f64le, filespace_pos,    dsetid_pos,    error)
        call h5dcreate_f(file_id, dsetname_index,  h5t_ieee_f64le, filespace_index,  dsetid_index,  error)
        call h5sclose_f(filespace_id,     error)
        call h5sclose_f(filespace_mass,   error)
        call h5sclose_f(filespace_charge, error)
        call h5sclose_f(filespace_vel,    error)
        call h5sclose_f(filespace_pos,    error)
        call h5sclose_f(filespace_index,  error)

        !> each process defines dataset in memory and writes it to the hyperslab in the file. 
        count_id(1)       = dim_id(1)
        count_id(2)       = count_part
        offset_id(1)      = 0
        offset_id(2)      = part_offset
        count_mass(1)     = dim_mass(1)
        count_mass(2)     = count_part
        offset_mass(1)    = 0
        offset_mass(2)    = part_offset
        count_charge(1)   = dim_charge(1)
        count_charge(2)   = count_part
        offset_charge(1)  = 0
        offset_charge(2)  = part_offset
        count_vel(1)      = dim_vel(1)
        count_vel(2)      = count_part
        offset_vel(1)     = 0
        offset_vel(2)     = part_offset
        count_pos(1)      = dim_pos(1)
        count_pos(2)      = count_part
        offset_pos(1)     = 0
        offset_pos(2)     = part_offset
        count_index(1)    = dim_index(1)
        count_index(2)    = count_part
        offset_index(1)   = 0
        offset_index(2)   = part_offset

        call h5screate_simple_f(rank_id,     count_id,     memspace_id,     error) 
        call h5screate_simple_f(rank_mass,   count_mass,   memspace_mass,   error) 
        call h5screate_simple_f(rank_charge, count_charge, memspace_charge, error) 
        call h5screate_simple_f(rank_vel,    count_vel,    memspace_vel,    error) 
        call h5screate_simple_f(rank_pos,    count_pos,    memspace_pos,    error) 
        call h5screate_simple_f(rank_index,  count_index,  memspace_index,  error) 

        !> select hyperslab in the file.
        call h5dget_space_f(dsetid_id,     filespace_id,     error)
        call h5dget_space_f(dsetid_mass,   filespace_mass,   error)
        call h5dget_space_f(dsetid_charge, filespace_charge, error)
        call h5dget_space_f(dsetid_vel,    filespace_vel,    error)
        call h5dget_space_f(dsetid_pos,    filespace_pos,    error)
        call h5dget_space_f(dsetid_index,  filespace_index,  error)
        call h5sselect_hyperslab_f(filespace_id,     h5s_select_set_f, offset_id,     count_id,     error)
        call h5sselect_hyperslab_f(filespace_mass,   h5s_select_set_f, offset_mass,   count_mass,   error)
        call h5sselect_hyperslab_f(filespace_charge, h5s_select_set_f, offset_charge, count_charge, error)
        call h5sselect_hyperslab_f(filespace_vel,    h5s_select_set_f, offset_vel,    count_vel,    error)
        call h5sselect_hyperslab_f(filespace_pos,    h5s_select_set_f, offset_pos,    count_pos,    error)
        call h5sselect_hyperslab_f(filespace_index,  h5s_select_set_f, offset_index,  count_index,  error)

        !> allocate data
        allocate(data_id(count_id(1), count_id(2)))
        allocate(data_mass(count_mass(1), count_mass(2)))
        allocate(data_charge(count_charge(1), count_charge(2)))
        allocate(data_vel(count_vel(1), count_vel(2)))
        allocate(data_pos(count_pos(1), count_pos(2)))
        allocate(data_index(count_index(1), count_index(2)))

        !> initilize data
        do i = 1, count_part
            data_id(1, i) = local_part_list(i) % id
            data_mass(1, i) = local_part_list(i) % mass + 1
            data_charge(1, i) = local_part_list(i) % charge
            data_vel(1, i) = local_part_list(i) % vel_x
            data_vel(2, i) = local_part_list(i) % vel_y
            data_vel(3, i) = local_part_list(i) % vel_z
            data_pos(1, i) = local_part_list(i) % pos_x
            data_pos(2, i) = local_part_list(i) % pos_y
            data_pos(3, i) = local_part_list(i) % pos_z
            data_index(1, i) = local_part_list(i) % global_cell_index_x
            data_index(2, i) = local_part_list(i) % global_cell_index_y
            data_index(3, i) = local_part_list(i) % global_cell_index_z
        end do

        ! !> create property list for collective dataset write
        call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error) 
        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, error)
     
        ! !> write the dataset collectively. 
        call h5dwrite_f(dsetid_id, h5t_native_double, data_id, dim_id, error, &
                        file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
        call h5dwrite_f(dsetid_mass, h5t_native_double, data_mass, dim_mass, error, &
                        file_space_id = filespace_mass, mem_space_id = memspace_mass, xfer_prp = plist_id)
        call h5dwrite_f(dsetid_charge, h5t_native_double, data_charge, dim_charge, error, &
                        file_space_id = filespace_charge, mem_space_id = memspace_charge, xfer_prp = plist_id)
        call h5dwrite_f(dsetid_vel, h5t_native_double, data_vel, dim_vel, error, &
                        file_space_id = filespace_vel, mem_space_id = memspace_vel, xfer_prp = plist_id)
        call h5dwrite_f(dsetid_pos, h5t_native_double, data_pos, dim_pos, error, &
                        file_space_id = filespace_pos, mem_space_id = memspace_pos, xfer_prp = plist_id)
        call h5dwrite_f(dsetid_index, h5t_native_double, data_index, dim_index, error, &
                        file_space_id = filespace_index, mem_space_id = memspace_index, xfer_prp = plist_id)

        deallocate(data_id)
        deallocate(data_mass)

        !> close dataspaces.
        call h5sclose_f(filespace_id,     error)
        call h5sclose_f(filespace_mass,   error)
        call h5sclose_f(filespace_charge, error)
        call h5sclose_f(filespace_vel,    error)
        call h5sclose_f(filespace_pos,    error)
        call h5sclose_f(filespace_index,  error)
        call h5sclose_f(memspace_id,      error)
        call h5sclose_f(memspace_mass,    error)
        call h5sclose_f(memspace_charge,  error)
        call h5sclose_f(memspace_vel,     error)
        call h5sclose_f(memspace_pos,     error)
        call h5sclose_f(memspace_index,   error)

        !> close the dataset and property list.
        call h5dclose_f(dsetid_id,     error)
        call h5dclose_f(dsetid_mass,   error)
        call h5dclose_f(dsetid_charge, error)
        call h5dclose_f(dsetid_vel,    error)
        call h5dclose_f(dsetid_pos,    error)
        call h5dclose_f(dsetid_index,  error)
        call h5pclose_f(plist_id, error)

        !> close the file.
        call h5fclose_f(file_id, error)

        !> close fortran predefined datatypes.
        call h5close_f(error)

    end subroutine output


    !> calulcate offset of output hyperslab
    function calculate_offset(counts_part, id)
        implicit none
        integer, dimension(numprocs)  :: counts_part
        integer                       :: id 
        integer                       :: calculate_offset
        integer                       :: i

        calculate_offset = 0
        do i = 0, id-1
            calculate_offset = calculate_offset + counts_part(i+1)
        end do
    end function


end module diagnostics