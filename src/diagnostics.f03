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
        character(len=8)                                       :: default_dir='../data'
        character(len=15)                                      :: filename
        !> dataset name
        character(len=13), parameter                           :: dsetname = "particle_info" 
        !> file identifier 
        integer(hid_t)                                         :: file_id       
        !> dataset identifier 
        integer(hid_t)                                         :: dset_id       
        !> dataspace identifier in file 
        integer(hid_t)                                         :: filespace     
        !> dataspace identifier in memory
        integer(hid_t)                                         :: memspace      
        !> property list identifier 
        integer(hid_t)                                         :: plist_id      

        !> dataset dimensions.
        integer(hsize_t), dimension(2)                         :: dimsf = (/5,8/) 
        integer(hsize_t), dimension(2)                         :: dimsfi = (/5,8/)

        integer(hsize_t), dimension(2)                         :: count  
        integer(hssize_t), dimension(2)                        :: offset 
        !> data to write
        double precision, allocatable                          :: data (:,:)  
        !> dataset rank 
        integer                                                :: rank = 2 

        integer                                                :: num_part(4)=(/ 1, 2, 2, 3 /)
        integer                                                :: offset_part(4)=(/ 0, 1, 3, 5 /)
        integer                                                :: i

        integer                                                :: error, error_n

        !> mpi related
        integer                                                :: comm, info
        integer                                                :: mpi_size, mpi_rank

        !> local particle
        integer, dimension(numprocs_x, numprocs_y, numprocs_z) :: locals_part_num
        integer                                                :: global_part_num

        call mpi_comm_size(cart_comm_3d, numprocs, ierr)
        call mpi_comm_rank(cart_comm_3d, my_id, ierr)

        !> output file name
        write(file_number, fmt) num
        filename = default_dir // trim(file_number) // '.h5'

        locals_part_num = count_local_particle(ierr)
        global_part_num = count_global_particle(ierr)


    end subroutine output

end module diagnostics