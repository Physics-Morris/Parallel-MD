module shared_data
    implicit none

    !> program control related shared data
    character(len=*), parameter :: VERSION = '1.0'
    logical                     :: green_light

    !> simulation domain related shared data
    integer                     :: sim_dimension
    double precision            :: x_min, x_max, y_min, y_max, z_min, z_max

    !> particle related shared data
    integer                     :: total_particles
    double precision            :: particle_mass, particle_charge
    character(len=60)           :: particle_distribution
    character(len=60)           :: velocity_distribution

    integer                     :: number_snapshots

    !> mpi related shared_data
    integer                     :: my_id, ierr, numprocs

end module shared_data
