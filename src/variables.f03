module variables
    implicit none

    character(len=*), parameter :: VERSION = '1.0'

    logical                     :: green_light

    integer     :: sim_dimension

    integer     :: total_particles

    integer     :: number_snapshots

    !> mpi related variables
    integer     :: my_id, ierr

    !> define particle type
    type particle
        integer          :: id, species_id, proc_id
        double precision :: mass, charge
        double precision :: pos_x, pos_y, pos_z
        double precision :: vel_x, vel_y, vel_z
    end type particle


end module variables
