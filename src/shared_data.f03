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

    !> rotate target
    logical                     :: rotate_target=.false.
    character(len=60)           :: rotate_sequence='xyz'
    character(len=60)           :: rotate_unit='degree'
    double precision            :: rotate_x=0.d0
    double precision            :: rotate_y=0.d0
    double precision            :: rotate_z=0.d0

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
    logical                     :: load_balance=.false.
    integer                     :: load_balance_num_step=-1
    double precision            :: load_balance_extent
    double precision            :: load_balance_threshold=-1.d0
    double precision            :: current_threshold
    double precision            :: max_speedup
    !> auxiliary cells
    integer, allocatable        :: auxi_cell(:, :, :, :)
    !> new auxi_cell after load balance
    integer, allocatable        :: auxi_cell_new(:, :, :, :)
    !> if load balance extent=0, use 10 auxi per processor
    integer                     :: num_auxi_per_procs=10
    integer                     :: max_num_auxi_per_procs=100
    integer                     :: auxi_num_x, auxi_num_y, auxi_num_z
    double precision            :: auxi_cell_wx, auxi_cell_wy, auxi_cell_wz

    !> z location of the dlb cut
    integer, allocatable        :: z_slice(:)
    !> y location of the dlb cut for particular (z) procs
    integer, allocatable        :: y_slice(:, :)
    !> x location of the dlb cut for particular (z, y) procs
    integer, allocatable        :: x_slice(:, :, :)

    
    !> initial distribution options parameter
    double precision            :: slab_xmin, slab_ymin, slab_zmin
    double precision            :: slab_xmax, slab_ymax, slab_zmax
    double precision            :: slab1_xmin, slab1_ymin, slab1_zmin
    double precision            :: slab1_xmax, slab1_ymax, slab1_zmax
    double precision            :: slab2_xmin, slab2_ymin, slab2_zmin
    double precision            :: slab2_xmax, slab2_ymax, slab2_zmax
    double precision            :: sphere_center_x, sphere_center_y, sphere_center_z
    double precision            :: sphere_radius
    double precision            :: FWHM_x, FWHM_y, FWHM_z
    double precision            :: slab_center_x, slab_center_y, slab_center_z


end module shared_data
