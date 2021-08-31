module shared_data
    use mpi
    implicit none

    !> program control related shared data
    character(len=*), parameter :: VERSION = 'v.0.1'
    logical                     :: green_light
    double precision            :: sim_start_time, sim_end_time

    !> simulation domain related shared data
    integer                     :: sim_dimension
    double precision            :: x_min, x_max, y_min, y_max, z_min, z_max
    !> local simulation domain
    double precision            :: x_min_local, x_max_local
    double precision            :: y_min_local, y_max_local
    double precision            :: z_min_local, z_max_local
    !> boundary condition
    character(len=20)           :: part_boundary_x, part_boundary_y, part_boundary_z

    integer                     :: step=0

    !> particle related shared data
    integer                     :: total_particles
    double precision            :: particle_mass, particle_charge
    character(len=60)           :: particle_distribution
    character(len=60)           :: velocity_distribution
    double precision            :: particle_temp_x, particle_temp_y, particle_temp_z
    integer                     :: local_particles

    integer                     :: number_snapshots

    !> mpi related shared_data
    double precision            :: task_start_time, task_end_time
    integer                     :: my_id, ierr, numprocs, cart_comm_3d
    integer                     :: particle_struc
    !> number of processor in each direction
    integer                     :: numprocs_x, numprocs_y, numprocs_z
    !> custome initial number of processor in each direction
    integer                     :: init_numprocs_x, init_numprocs_y, init_numprocs_z
    !> load balance option
    logical                     :: load_balance
    integer                     :: load_balance_num_step
    !> auxiliary cells
    integer, allocatable        :: auxi_cell(:, :, :, :)
    integer                     :: num_auxi_per_procs=100
    integer                     :: auxi_num_x, auxi_num_y, auxi_num_z
    double precision            :: auxi_cell_wx, auxi_cell_wy, auxi_cell_wz

end module shared_data
