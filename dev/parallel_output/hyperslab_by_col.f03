
! number of processes is assumed to be 1 or multiples of 2 (1,2,4,6,8)
!

     program dataset_by_col

     use hdf5 ! this module contains all necessary modules 
     use mpi
        
     implicit none

     ! include 'mpif.h'
     character(len=10), parameter :: filename = "sds_col.h5"  ! file name
     character(len=13), parameter :: dsetname = "particle_info" ! dataset name
     character(len=8), parameter :: dsetname2 = "velocity" ! dataset name

     integer(hid_t) :: file_id       ! file identifier 
     integer(hid_t) :: dset_id       ! dataset identifier 
     integer(hid_t) :: filespace     ! dataspace identifier in file 
     integer(hid_t) :: memspace      ! dataspace identifier in memory
     integer(hid_t) :: plist_id      ! property list identifier 

     integer(hid_t) :: dset_id2       ! dataset identifier 
     integer(hid_t) :: filespace2     ! dataspace identifier in file 
     integer(hid_t) :: memspace2      ! dataspace identifier in memory

     integer(hsize_t), dimension(2) :: dimsf = (/5,8/) ! dataset dimensions.
!     integer, dimension(7) :: dimsfi = (/5,8,0,0,0,0,0/)
     integer(hsize_t), dimension(2) :: dimsfi = (/5,8/)

     integer(hsize_t), dimension(2) :: count  
     integer(hssize_t), dimension(2) :: offset 
     double precision, allocatable :: data (:,:)  ! data to write
     integer :: rank = 2 ! dataset rank 

     integer   :: num_part(4)=(/ 1, 2, 2, 3 /)
     integer   :: offset_part(4)=(/ 0, 1, 3, 5 /)
     integer   :: i

     integer :: error, error_n  ! error flags
     !
     ! mpi definitions and calls.
     !
     integer :: mpierror       ! mpi error flag
     integer :: comm, info
     integer :: mpi_size, mpi_rank

     comm = mpi_comm_world
     info = mpi_info_null
     call mpi_init(mpierror)
     call mpi_comm_size(comm, mpi_size, mpierror)
     call mpi_comm_rank(comm, mpi_rank, mpierror) 
     !
     ! initialize fortran predefined datatypes
     !
     call h5open_f(error) 

     ! 
     ! setup file access property list with parallel i/o access.
     !
     call h5pcreate_f(h5p_file_access_f, plist_id, error)
     call h5pset_fapl_mpio_f(plist_id, comm, info, error)

     !
     ! create the file collectively.
     ! 
     call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error, access_prp = plist_id)
     call h5pclose_f(plist_id, error)
     !
     ! create the data space for the  dataset. 
     !
     call h5screate_simple_f(rank, dimsf, filespace, error)
     call h5screate_simple_f(rank, dimsf, filespace2, error)

     !
     ! create the dataset with default properties.
     !
     call h5dcreate_f(file_id, dsetname, h5t_ieee_f64le, filespace, &
                      dset_id, error)
     call h5sclose_f(filespace, error)
     call h5dcreate_f(file_id, dsetname2, h5t_ieee_f64le, filespace2, &
                      dset_id2, error)
     call h5sclose_f(filespace2, error)
     !
     ! each process defines dataset in memory and writes it to the hyperslab
     ! in the file. 
     !
     count(1) = dimsf(1)
     ! count(2) = dimsf(2)/mpi_size 
     count(2) = num_part(mpi_rank+1)
     offset(1) = 0
     ! offset(2) = mpi_rank * count(2) 
     offset(2) = offset_part(mpi_rank+1)
     call h5screate_simple_f(rank, count, memspace, error) 
     call h5screate_simple_f(rank, count, memspace2, error) 
     ! 
     ! select hyperslab in the file.
     !
     call h5dget_space_f(dset_id, filespace, error)
     call h5sselect_hyperslab_f (filespace, h5s_select_set_f, offset, count, error)
     call h5dget_space_f(dset_id2, filespace2, error)
     call h5sselect_hyperslab_f (filespace2, h5s_select_set_f, offset, count, error)
     ! 
     ! initialize data buffer with trivial data.
     !
     allocate ( data(count(1),count(2)))
     do i = 1, count(1)
          data(i, :) = dble(mpi_rank) + 10.d0 + dble(i)/10.d0
     end do
     !
     ! create property list for collective dataset write
     !
     call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error) 
     call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, error)
     
     !
     ! write the dataset collectively. 
     !
     call h5dwrite_f(dset_id, h5t_native_double, data, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
     call h5dwrite_f(dset_id2, h5t_native_double, data, dimsfi, error, &
                     file_space_id = filespace2, mem_space_id = memspace2, xfer_prp = plist_id)
     !
     ! write the dataset independently. 
     !
!    call h5dwrite_f(dset_id, h5t_native_integer, data, dimsfi, error, &
!                     file_space_id = filespace, mem_space_id = memspace)
     !
     ! deallocate data buffer.
     !
     deallocate(data)

     !
     ! close dataspaces.
     !
     call h5sclose_f(filespace, error)
     call h5sclose_f(memspace, error)
     call h5sclose_f(filespace2, error)
     call h5sclose_f(memspace2, error)

     !
     ! close the dataset and property list.
     !
     call h5dclose_f(dset_id, error)
     call h5dclose_f(dset_id2, error)
     call h5pclose_f(plist_id, error)

     !
     ! close the file.
     !
     call h5fclose_f(file_id, error)

     !
     ! close fortran predefined datatypes.
     !
     call h5close_f(error)

     call mpi_finalize(mpierror)

     end program dataset_by_col