module particle
    use constants
    use shared_data
    use math
    implicit none

    !> derived particle type
    type, public :: particle_type
        integer          :: id
        integer          :: global_cell_index_x
        integer          :: global_cell_index_y
        integer          :: global_cell_index_z
        double precision :: mass, charge
        double precision :: pos_x, pos_y, pos_z
        double precision :: vel_x, vel_y, vel_z

        contains

        !> new features in Fortran 2003, allow us to have type bound procedure
        procedure :: kinetic_energy => calculate_kinetic_energy
        procedure :: gamma => calculate_gamma
        procedure :: speed => calculate_speed
        procedure :: rest_energy => calculate_rest_energy
        procedure :: momentum => calculate_momentum
        procedure :: initialize => initialize_particle
    end type particle_type


    !> the list that store global partarticle and local ones
    type(particle_type), dimension(:), allocatable :: global_part_list
    type(particle_type), dimension(:), allocatable :: local_part_list


    contains

    !> caluclate speed of particle
    function calculate_speed(this) result(speed)
        implicit none
        class(particle_type), intent(in) :: this
        double precision :: speed
            speed = dsqrt(this%vel_x**2 + this%vel_y**2 + this%vel_z**2)
    end function calculate_speed


    !> calculate gamma coefficient of the particle
    function calculate_gamma(this) result(gamma)
        implicit none
        class(particle_type), intent(in) :: this
        double precision  ::  gamma
        gamma = 1.d0 / dsqrt(1.d0 - calculate_speed(this)**2/c**2)
    end function calculate_gamma


    !> calculate kinetic energy of the particle
    function calculate_kinetic_energy(this) result(KE)
        implicit none
        class(particle_type), intent(in) :: this
        double precision  ::  KE
        !> relativistic kinetic energy
        KE = (calculate_gamma(this) - 1.d0) * this%mass * c**2 
    end function calculate_kinetic_energy


    !> calculate rest energy of the particle
    function calculate_rest_energy(this) result(rest_energy)
        implicit none
        class(particle_type), intent(in) :: this
        double precision  ::  rest_energy
        !> rest energy
        rest_energy = this%mass * c**2
    end function calculate_rest_energy


    !> calculate relativistic momentum of the particle
    function calculate_momentum(this) result(p)
        implicit none
        class(particle_type), intent(in) :: this
        double precision  ::  p
        p = calculate_gamma(this) * this%mass * calculate_speed(this)
    end function calculate_momentum


    !> initialize particle derived type
    subroutine initialize_particle(this)
        implicit none
        class(particle_type), intent(inout) :: this
        this % id = 0
        this % global_cell_index_x = 0
        this % global_cell_index_y = 0
        this % global_cell_index_z = 0
        this % mass = 0.d0
        this % charge = 0.d0
        this % pos_x = 0.d0
        this % pos_y = 0.d0
        this % pos_z = 0.d0
        this % vel_x = 0.d0
        this % vel_y = 0.d0
        this % vel_z = 0.d0
    end subroutine initialize_particle


    !> load particle into simulation space globally
    subroutine load_particles_globally(ierr)
        implicit none
        integer, intent(inout) :: ierr
        integer                :: i
        
        !> first allocate global particle list
        allocate(global_part_list(total_particles), stat=ierr)

        !> initialize all particles
        do i = 1, total_particles
            call global_part_list(i) % initialize
            !> asign id
            global_part_list(i) % id = i
        end do

        !> first load particle position distribution accroding to
        !> density function
        if (particle_distribution == 'uniform_in_domain') then
            do i = 1, total_particles
                global_part_list(i) % pos_x = uniform_distribution(x_min, x_max)
                global_part_list(i) % pos_y = uniform_distribution(y_min, y_max)
                global_part_list(i) % pos_z = uniform_distribution(z_min, z_max)
            end do

        else if (particle_distribution == 'custom') then
            !>
            ! use custom fucntion locate in user_custom directory
            ! (add in the future)
            !>
        end if

        !> load particle velocity distrituion
        if (velocity_distribution == 'cold_xyz') then
            do i = 1, total_particles
                global_part_list(i) % vel_x = 0.d0
                global_part_list(i) % vel_y = 0.d0
                global_part_list(i) % vel_z = 0.d0
            end do
        end if
    end subroutine load_particles_globally


    !> unload particle globally
    subroutine unload_particles_globally(ierr)
        implicit none
        integer, intent(inout) :: ierr
        
        deallocate(global_part_list, stat=ierr)
    end subroutine unload_particles_globally


    !> unload particle locally
    subroutine unload_particles_locally(ierr)
        implicit none
        integer, intent(inout) :: ierr
        
        deallocate(local_part_list, stat=ierr)
    end subroutine unload_particles_locally


    !> map particle in global cell
    function map_particle_to_global_cell(pos_x, pos_y, pos_z)
        implicit none
        integer, dimension(3)        :: map_particle_to_global_cell
        double precision, intent(in) :: pos_x, pos_y, pos_z
        double precision             :: wx, wy, wz

        wx = (x_max - x_min) / dble(numprocs_x)
        wy = (y_max - y_min) / dble(numprocs_y)
        wz = (z_max - z_min) / dble(numprocs_z)

        map_particle_to_global_cell(1) = floor((pos_x - x_min) / wx) + 1
        map_particle_to_global_cell(2) = floor((pos_y - y_min) / wy) + 1
        map_particle_to_global_cell(3) = floor((pos_z - z_min) / wz) + 1

    end function map_particle_to_global_cell



end module particle