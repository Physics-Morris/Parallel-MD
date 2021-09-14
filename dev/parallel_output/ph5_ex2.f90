PROGRAM main

  USE HDF5
  USE ISO_C_BINDING

  IMPLICIT NONE
  INCLUDE 'mpif.h'
  CHARACTER(LEN=17), PARAMETER :: filename  = "h5ex_t_int_F03.h5"
  CHARACTER(LEN=2) , PARAMETER :: dataset   = "DS"
  INTEGER          , PARAMETER :: dim0      = 4
  INTEGER          , PARAMETER :: dim1      = 7
  INTEGER          , PARAMETER :: ng        = 64
  INTEGER          , PARAMETER :: nd        = 16

  INTEGER(HID_T)  :: file, space, dset, gid, gid1, plist_id ! Handles

  INTEGER(hsize_t),   DIMENSION(1:2) :: dims = (/dim0, dim1/)
  INTEGER, DIMENSION(1:dim0, 1:dim1), TARGET :: wdata ! Write buffer
  INTEGER, DIMENSION(:,:), ALLOCATABLE, TARGET :: rdata ! Read buffer
  INTEGER(HSIZE_T), DIMENSION(1:2) :: maxdims
  INTEGER :: i, j, jng, jnd, k
  TYPE(C_PTR) :: f_ptr
  CHARACTER*3 :: ichr3
  INTEGER :: ierr
  INTEGER :: mpierror
  INTEGER :: mpi_size, mpi_rank
  
  !
  ! Initialize FORTRAN interface.
  !

  CALL MPI_INIT(mpierror)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, mpierror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, mpierror) 
  
  CALL h5open_f(ierr)
  !
  ! Initialize DATA.
  !
  DO i = 1, dim0
     DO j = 1, dim1
        wdata(i,j) = (i-1) * (j-1) - (j-1)
     ENDDO
  ENDDO

  IF( mpi_rank .EQ. 0) THEN
     
     !
     ! Create a new file using the default properties.
     !
     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, ierr)
     !
     ! Create dataspace.  Setting maximum size to be the current size.
     !
     CALL h5screate_simple_f(2, dims, space, ierr)
     !
     ! Create the dataset and write the integer data to it.  In this
     ! example we will save the data as 64 bit big endian integers,
     ! regardless of the native integer type.  The HDF5 library
     ! automatically converts between different integer types.
     !

     f_ptr = C_LOC(wdata(1,1))
     DO i = 1, ng
        WRITE(ichr3, '(I3.3)') i
        CALL h5gcreate_f(file,"grp_"//ichr3, gid, ierr)
        CALL h5gcreate_f(gid,"grp_"//ichr3, gid1, ierr)
        DO j = 1, nd
           WRITE(ichr3, '(I3.3)') j
           CALL h5dcreate_f(gid1, TRIM(dataset)//ichr3, H5T_STD_I64BE, space, dset, ierr)
           CALL h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, ierr)
           CALL h5dclose_f(dset , ierr)
        ENDDO
        
        CALL h5gclose_f(gid,ierr)
        CALL h5gclose_f(gid1,ierr)
     ENDDO
     
     !
     ! Close and release resources.
     !
     CALL h5sclose_f(space, ierr)
     CALL h5fclose_f(file , ierr)
  ENDIF
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
  CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)

  !
  ! Open file and dataset.
  !
  CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file, ierr, plist_id)
  CALL h5pclose_f(plist_id, ierr)

  jng = ng/mpi_size
  jnd = nd/mpi_size

  DO i = jng*mpi_rank+1, jng*(1+mpi_rank)
     WRITE(ichr3, '(I3.3)') i
     CALL h5gopen_f(file,"grp_"//ichr3, gid, ierr)
     CALL h5gopen_f(gid,"grp_"//ichr3, gid1, ierr)
     DO j = 1, nd
        WRITE(ichr3, '(I3.3)') j
        CALL h5dopen_f(gid1, dataset//ichr3, dset, ierr)    
        !
        ! Get dataspace and allocate memory for read buffer.
        !
        CALL h5dget_space_f(dset,space, ierr)
        CALL h5sget_simple_extent_dims_f(space, dims, maxdims, ierr)

        ALLOCATE(rdata(1:dims(1),1:dims(2)))
        !
        ! Read the data.
        !
        f_ptr = C_LOC(rdata(1,1))
        CALL h5dread_f( dset, H5T_NATIVE_INTEGER, f_ptr, ierr)
        !
        ! Output the data to the screen.
        !
#if 0
        IF(mpi_rank.EQ.2)THEN
           WRITE(*, '(A,":")') dataset//ichr3
           DO k=1, dims(1)
              WRITE(*,'(" [")', ADVANCE='NO')
              WRITE(*,'(80i3)', ADVANCE='NO') rdata(k,1:dims(2))
              WRITE(*,'(" ]")')
           ENDDO
        ENDIF
#endif
  !
  ! Close and release resources.
  !
        DEALLOCATE(rdata)
        CALL h5dclose_f(dset , ierr)
        CALL h5sclose_f(space, ierr)
     END DO
     CALL h5gclose_f(gid, ierr)
     CALL h5gclose_f(gid1, ierr)


  ENDDO

  CALL h5fclose_f(file , ierr)
  CALL h5close_f(ierr)

  CALL MPI_FINALIZE(mpierror)

END PROGRAM main
