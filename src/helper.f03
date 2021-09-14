module helper
    use shared_data
    use mpi
    use constants
    use math
    use omp_lib
    use, intrinsic :: iso_c_binding, only: c_int, c_int32_t
    implicit none


    namelist / basics_block / &
    sim_dimension, x_min, x_max, y_min, y_max, z_min, z_max, &
    init_numprocs_x, init_numprocs_y, init_numprocs_z, load_balance, &
    load_balance_num_step, load_balance_extent, load_balance_threshold

    namelist / particles_block / &
    total_particles, particle_mass, particle_charge, particle_distribution, &
    slab_xmin, slab_xmax, slab_ymin, slab_ymax, slab_zmin, slab_zmax, &
    slab1_xmin, slab1_xmax, slab1_ymin, slab1_ymax, slab1_zmin, slab1_zmax, &
    slab2_xmin, slab2_xmax, slab2_ymin, slab2_ymax, slab2_zmin, slab2_zmax, &
    sphere_center_x, sphere_center_y, sphere_center_z, sphere_radius, &
    FWHM_x, FWHM_y, FWHM_z, slab_center_x, slab_center_y, slab_center_z, &
    velocity_distribution, particle_temp_x, particle_temp_y, particle_temp_z, &
    part_boundary_x, part_boundary_y, part_boundary_z, rotate_target, rotate_sequence, &
    rotate_unit, rotate_x, rotate_y, rotate_z

    namelist / output_block / &
    number_snapshots

    interface
        ! int usleep(useconds_t useconds)
        function c_usleep(useconds) bind(c, name='usleep')
            import :: c_int, c_int32_t
            integer(kind=c_int32_t), value :: useconds
            integer(kind=c_int)            :: c_usleep
        end function c_usleep
    end interface

    contains


        subroutine get_input_file(loc, ierr)
            implicit none
            character(len=*)            :: loc
            logical                     :: lexist
            integer                     :: ierr

            !> inquire input file
            call print_empty_line(1)
            call set_term_color(term_default_colour)
            call print_execute_task_name('  Inquire input file ... ', task_start_time)
            inquire (file=trim(loc), exist=lexist)
            if (lexist) then
                call print_task_time(task_start_time)
                call successful
                call set_term_color(term_default_colour)
                call print_empty_line(1)
            else
                call print_task_time(task_start_time)
                call failed
                call set_term_color(term_default_colour)
                call print_empty_line(1)
                write(*, *) ' Input file does not exist, stopping the program ... '
                call print_empty_line(1)
                stop
            end if

            !> read input file
            write(*, '(A, A, A)', advance='no') '  Reading input file (', trim(loc)
            call print_execute_task_name(') ... ', task_start_time)
            open(10, file=trim(loc), status='old', action='read', iostat=ierr)
            if (ierr .eq. 0) then
                read(10, nml=basics_block)
                read(10, nml=particles_block)
                read(10, nml=output_block)
                call print_task_time(task_start_time)
                call successful
                close(10)
            else
                call print_task_time(task_start_time)
                call failed
                write(*, '(A, I5)') '  error message number: ', ierr
                stop
            end if
            call print_simulation_parameter('full')
            call sleep(1)
            
        end subroutine get_input_file


        subroutine print_simulation_parameter(status)
            implicit none
            character(len=*)             :: status

            if (status == 'full') then
                call print_empty_line(2)
                call print_divider('-', 55)

                call set_term_color(term_default_colour)

                call check_input_parameter('simulation_dimension')
                call check_input_parameter('init_procs_layout')
                call check_input_parameter('load_balance')
                call check_input_parameter('load_balance_step')
                call check_input_parameter('load_balance_extent')
                call check_input_parameter('load_balance_threshold')
                call check_input_parameter('x_min')
                call check_input_parameter('x_max')
                call check_input_parameter('y_min')
                call check_input_parameter('y_max')
                call check_input_parameter('z_min')
                call check_input_parameter('z_max')
                call check_input_parameter('total_particles')
                call check_input_parameter('particle_mass')
                call check_input_parameter('particle_charge')
                call check_input_parameter('particle_distribution')
                call check_input_parameter('rotate_target')
                call check_input_parameter('particle_velocity_distribution')
                call check_input_parameter('temp_x')
                call check_input_parameter('temp_y')
                call check_input_parameter('temp_z')
                call check_input_parameter('part_boundary_x')
                call check_input_parameter('part_boundary_y')
                call check_input_parameter('part_boundary_z')
                call check_input_parameter('number_output_snapshots')

                call set_term_color(term_default_colour)
                call print_divider('-', 55)
                call print_empty_line(2)
            end if
        end subroutine print_simulation_parameter


        subroutine get_cmd_arg
            implicit none
            character(len=32)           :: arg
            integer                     :: i, ierr

            !> must specify input file
            if (command_argument_count() .eq. 0) then
                CALL set_term_color(3)
                call print_empty_line(1)
                write(*, *) ' Need to specify input file'
                call print_empty_line(1)
                CALL set_term_color(12)
                call print_help()
                stop
            end if

            do i = 1, command_argument_count()
                call get_command_argument(i, arg)
                select case (arg)
                    case ('-v', '--version')
                        call print_empty_line(2)
                        print '(2a)', '             version ', VERSION
                        call MD_logo
                        stop
        
                    case ('-h', '--help')
                        call print_help()
                        stop

                    case ('-i', '--input')
                        call get_command_argument(i+1, arg)
                        call get_input_file(arg, ierr)
                        exit

                    case ('-t', '--time')
                        call print_time
                        stop

                    case ('-d', '--default')
                        call print_empty_line(1)
                        write(*, *) ' Useing default input file... '
                        call get_input_file('../inp/default.input', ierr)
                        exit

                    case default
                        print '(2a, /)', ' unrecognised command-line option: ', arg
                        call print_help()
                        stop
                end select 
            end do
        end subroutine get_cmd_arg

        subroutine print_time
            character(len=8) :: date
            character(len=10) :: time
            character(len=5) :: zone

            !> Print the date and time
            call set_term_color(5)
            call print_empty_line(1)
            call date_and_time(DATE=date, TIME=time, ZONE=zone)
            write (*, '(a, a, "-", a, "-", a, a, a, ":", a)') ' Current date and time: ', &
            date(1:4), date(5:6), date(7:8), '  ', time(1:2), time(3:4)
            call print_empty_line(1)
            call set_term_color(12)
        end subroutine print_time

        subroutine print_help
            implicit none
            call set_term_color(9)
            write(*, '(a, /)')   '                          command-line options'
            write(*, *)          ' ======================================================================='
            write(*, '(a)')      '  -i, --input            , entering input file location'
            write(*, '(a)')      '  -d, --default          , use default input file for testing'
            write(*, '(a)')      '  -v, --version          , check out the version of the program'
            write(*, '(a)')      '  -h, --help             , see hlep for possible command line arguments'
            write(*, '(a)')      '  -t, --time             , print current date and time'
            write(*, *)          ' ======================================================================='
            call set_term_color(12)
            call print_empty_line(1)
        end subroutine print_help


        subroutine welcome
            implicit none
            call MD_logo
            call set_term_color(term_green)
            call print_empty_line(1)
            CALL set_term_color(term_default_colour)
            call set_term_color(term_reset_attributes)
        end subroutine welcome

        
        subroutine MD_logo
            implicit none
            integer :: i
            call print_empty_line(1)
            CALL set_term_color(term_green)
            write(*, '(A)', advance='no') '  All greens, ready to go '
            do i = 1, 5
                write(*, '(A)', advance='no') '.'
                call sleep(1)
            end do
            CALL set_term_color(term_default_colour)
            call print_empty_line(1)
            call set_term_color(term_bold)
            call set_term_color(term_yellow)
            call print_empty_line(4)
            call set_term_color(3)
            write(*,'(A)') '            ###     ###     ######## '
            call set_term_color(5)
            write(*,'(A)') '           ####   ####     ######### '
            call set_term_color(2)
            write(*,'(A)') '          -----------     --      -- '
            call set_term_color(4)
            write(*,'(A)') '         ### ### ###     ##       ## '
            call set_term_color(6)
            write(*,'(A)') '        ###     ###     ##       ##  '
            call set_term_color(8)
            write(*,'(A)') '       ###     ###     ##       ##   '
            call set_term_color(7)
            write(*,'(A)') '      ###     ###     ##########     '
            call set_term_color(1)
            write(*,'(A)') '     ###     ###     #########       '
            call print_empty_line(1)
            call set_term_color(12)
        end subroutine MD_logo


        subroutine MICPIC_logo
            implicit none
            call print_empty_line(1)
            CALL set_term_color(term_green)
            write(*, '(A)') ' All greens, ready to go'
            CALL set_term_color(12)
            call sleep(3)
            CALL set_term_color(term_bold)
            CALL set_term_color(term_yellow)
            call print_empty_line(4)
            call set_term_color(3)
            write(*,'(A)') '            ###     ###                    #########     #########       #######'
            call set_term_color(1)
            write(*,'(A)') '           ####   ####                    ###########   #########     ######### '
            call set_term_color(2)
            call set_term_color(term_cyan)
            call set_term_color(term_yellow)
            write(*,'(A)') '          -----------   ---    -----     ----    ----      ---       -----      '
            call set_term_color(4)
            write(*,'(A)') '         ### ### ###   ####   ######    ###########       ###      ######       '
            call set_term_color(5)
            write(*,'(A)') '        ###     ###     ##   ###       #########         ###      #######       '
            call set_term_color(6)
            write(*,'(A)') '       ###     ###     ##   ###       #####             ###       ######        '
            call set_term_color(7)
            write(*,'(A)') '      ###     ###     ##    ###      #####          ##########     ########     '
            call set_term_color(8)
            write(*,'(A)') '     ###     ###     ###    #####   #####          ##########      #######      '
            call print_empty_line(1)
            call set_term_color(12)
        end subroutine MICPIC_logo


        subroutine print_empty_line(num)
            implicit none
            integer             :: i, num
            do i = 1, num
                write(*, *)
            end do
        end subroutine print_empty_line


        subroutine set_term_color(controlcode)
            implicit none
            integer, intent(in) :: controlcode
            write(*, '(a)', advance='no') achar(27) // trim(vt100_control(controlcode))
        end subroutine set_term_color


        subroutine successful
            implicit none
            integer :: rc
            rc = c_usleep(100*1000)
            call set_term_color(term_green)
            green_light = .True.
            write(*, *) 'Successful'
            call set_term_color(term_default_colour)
        end subroutine successful


        subroutine failed
            implicit none
            integer :: rc
            rc = c_usleep(100*1000)
            call set_term_color(term_red)
            green_light = .False.
            write(*, *) 'Failed'
            call set_term_color(term_default_colour)
        end subroutine failed


        subroutine respond_to_ierr(ierr)
            implicit none
            integer, intent(in) :: ierr
            call mpi_comm_rank(mpi_comm_world, my_id, ierr)
            if (my_id .eq. master_id) then
                if (ierr == 0) then
                    call successful
                else
                    call failed
                end if
                call print_empty_line(1)
            end if
        end subroutine respond_to_ierr


        subroutine print_execute_task_name(name, task_start_time)
            implicit none
            character(len=*) :: name
            double precision :: task_start_time
            call mpi_comm_rank(mpi_comm_world, my_id, ierr)
            if (my_id .eq. master_id) then
                task_start_time = mpi_wtime()
                call set_term_color(term_default_colour)
                write(*, '(A)', advance='no') name
            end if
        end subroutine print_execute_task_name
        

        subroutine print_task_time(task_start_time)
            implicit none
            double precision :: task_start_time, task_end_time
            call mpi_comm_rank(mpi_comm_world, my_id, ierr)
            if (my_id .eq. master_id) then
                task_end_time = mpi_wtime()
                call set_term_color(term_yellow)
                write(*, '(A, ES7.1, A)', advance='no') '(', task_end_time-task_start_time, ' sec)'
                call set_term_color(term_default_colour)
            end if
        end subroutine print_task_time
        

        subroutine sec_to_day_hour_min(day, hour, min, sec)
            implicit none
            integer, intent(out)   :: day, hour, min
            integer, intent(inout) :: sec
            day  = int(floor(dble(sec) / dble(86400)))
            hour = int(floor(dble(sec) / dble(3600)))
            min  = int((dble(sec) / dble(3600) - hour) * 60)
            sec  = int((((dble(sec) / dble(3600) - hour) * 60) - min) * 60)
        end subroutine


        subroutine print_ok_mark(ierr)
            implicit none
            integer :: ierr, rc
            if (ierr == 0) then
                rc = c_usleep(100*1000)
                call set_term_color(term_green)
                write(*, '(A3)') ' OK' 
                call set_term_color(term_default_colour)
                rc = c_usleep(100*1000)
            else
                rc = c_usleep(100*1000)
                call set_term_color(term_red)
                write(*, '(A2)') ' X' 
                call set_term_color(term_default_colour)
                rc = c_usleep(100*1000)
            end if
        end subroutine print_ok_mark


        subroutine check_input_parameter(name)
            implicit none
            character(len=*), intent(in)  :: name
            double precision              :: acpt_range=1.d-5
            integer                       :: default_length=25

            select case (name)

            case('simulation_dimension')
                !> currently only support 3d
                write(*, '(A30)', advance='no') '  Simulation dimension:       '
                if (sim_dimension == 3) then
                    write(*, '(I25)', advance='no') sim_dimension
                    call print_ok_mark(ierr=0)
                else 
                    write(*, '(I25)', advance='no') sim_dimension
                    call print_ok_mark(ierr=1)
                    stop
                end if

            case('init_procs_layout')
                !> sepcify proccsor layout must match numprocs
                call mpi_comm_size(mpi_comm_world, numprocs, ierr)
                write(*, '(A35)', advance='no') '  Initial proccesor layout:        '
                if (init_numprocs_x*init_numprocs_y*init_numprocs_z == numprocs) then
                    write(*, '(A6, A2, I2, A2, I2, A2, I2, A2)', advance='no') &
                    '      ', '[', init_numprocs_x, 'x', init_numprocs_y, 'x', init_numprocs_z, ']'
                    call print_ok_mark(ierr=0)
                else 
                    write(*, '(A2, I2, A2, I2, A2, I2, A2)', advance='no') &
                    '[', init_numprocs_x, 'x', init_numprocs_y, 'x', init_numprocs_z, ']'
                    call print_ok_mark(ierr=1)
                    stop
                end if

            case('load_balance')
                write(*, '(A30, L25)', advance='no')     '  Load balancing:             ', load_balance
                call print_ok_mark(ierr=0)

            case('load_balance_step')
                !> load balance step must greater than 0
                if ((load_balance .eqv. .true.).and.(load_balance_num_step >= 0)) then
                    write(*, '(A30)', advance='no') '  Load balancing every (step):'
                    write(*, '(I25)', advance='no') load_balance_num_step
                    call print_ok_mark(ierr=0)
                end if

            case('load_balance_extent')
                write(*, '(A30)', advance='no') '  Load balancing extent [0, 1]:'
                !> load balance step must greater than 0
                if ((load_balance_extent < 0).or.(load_balance_extent > 1)) then
                    write(*, '(F25.2)', advance='no') load_balance_extent
                    call print_ok_mark(ierr=1)
                    stop
                else 
                    num_auxi_per_procs = num_auxi_per_procs + &
                                         nint(max_num_auxi_per_procs*load_balance_extent)
                    write(*, '(F25.2)', advance='no') load_balance_extent
                    call print_ok_mark(ierr=0)
                end if

            case('load_balance_threshold')
                !> load balance step must greater than 0
                if ((load_balance .eqv. .true.).and.(load_balance_num_step < 0).and. &
                    (load_balance_threshold >= 0.d0).and.(load_balance_threshold <= 1.d0)) then
                    write(*, '(A30)', advance='no') '  Load threshold [0, 1]:      '
                    write(*, '(F25.2)', advance='no') load_balance_threshold
                    call print_ok_mark(ierr=0)
                else if ((load_balance .eqv. .true.).and.(load_balance_num_step < 0).and. &
                         (load_balance_threshold < 0.d0).or.(load_balance_threshold > 1.d0)) then
                    write(*, *) ' Needs to specify either load_balance_num_step or load_balance_threshold'
                    write(*, *) ' and choose right range of load_balance_num_step or load_balance_threshold'
                    stop
                else if ((load_balance .eqv. .true.).and.(load_balance_num_step >= 0).and. &
                         (load_balance_threshold >= 0.d0)) then
                    write(*, *) ' You can only choose either load_balance_num_step or load_balance_threshold'
                    stop
                end if

            case('x_min')
                write(*, '(A30)', advance='no') '  x_min of the system:        '
                if (x_min >= x_max) then
                    write(*, '(F25.2)', advance='no') x_min
                    call print_ok_mark(ierr=1)
                else 
                    write(*, '(F25.2)', advance='no') x_min
                    call print_ok_mark(ierr=0)
                end if

            case('x_max')
                write(*, '(A30)', advance='no') '  x_max of the system:        '
                if (x_min >= x_max) then
                    write(*, '(F25.2)', advance='no') x_max
                    call print_ok_mark(ierr=1)
                    stop
                else 
                    write(*, '(F25.2)', advance='no') x_max
                    call print_ok_mark(ierr=0)
                end if

            case('y_min')
                write(*, '(A30)', advance='no') '  y_min of the system:        '
                if (y_min >= y_max) then
                    write(*, '(F25.2)', advance='no') y_min
                    call print_ok_mark(ierr=1)
                else 
                    write(*, '(F25.2)', advance='no') y_min
                    call print_ok_mark(ierr=0)
                end if

            case('y_max')
                write(*, '(A30)', advance='no') '  y_max of the system:        '
                if (y_min >= y_max) then
                    write(*, '(F25.2)', advance='no') y_max
                    call print_ok_mark(ierr=1)
                    stop
                else 
                    write(*, '(F25.2)', advance='no') y_max
                    call print_ok_mark(ierr=0)
                end if

            case('z_min')
                write(*, '(A30)', advance='no') '  z_min of the system:        '
                if (z_min >= z_max) then
                    write(*, '(F25.2)', advance='no') z_min
                    call print_ok_mark(ierr=1)
                else 
                    write(*, '(F25.2)', advance='no') z_min
                    call print_ok_mark(ierr=0)
                end if

            case('z_max')
                write(*, '(A30)', advance='no') '  z_max of the system:        '
                if (z_min >= z_max) then
                    write(*, '(F25.2)', advance='no') z_max
                    call print_ok_mark(ierr=1)
                    stop
                else 
                    write(*, '(F25.2)', advance='no') z_max
                    call print_ok_mark(ierr=0)
                end if

            case('total_particles')
                write(*, '(A30, ES25.0)', advance='no')     '  Total number of particles:  ', dble(total_particles)
                if (total_particles <= 0) then
                    call print_ok_mark(ierr=1)
                    stop
                else 
                    call print_ok_mark(ierr=0)
                end if

            case('particle_mass')
                write(*, '(A30)', advance='no') '  Particle mass:              '
                !> if particle mass = -1 than it is electron
                if (abs(particle_mass+1.d0) <= acpt_range) then
                    particle_mass = m0
                    write(*, '(ES25.2)', advance='no') particle_mass
                    call print_ok_mark(ierr=0)
                !> if particle mass = 1 than it is proton
                else if (abs(particle_mass-1.d0) <= acpt_range) then
                    particle_mass = 1836.d0 * m0
                    write(*, '(ES25.2)', advance='no') particle_mass
                    call print_ok_mark(ierr=0)
                !> if particle mass < -1 than is it the number of times of
                !> the mass of electron
                else if ((particle_mass < 0.d0).and.(abs(particle_mass+1.d0) >= acpt_range)) then
                    particle_mass = abs(particle_mass) * m0
                    write(*, '(ES25.2)', advance='no') particle_mass
                    call print_ok_mark(ierr=0)
                else 
                    write(*, '(ES25.2)', advance='no') particle_mass
                    call print_ok_mark(ierr=0)
                end if

            case('particle_charge')
                write(*, '(A30)', advance='no') '  Particle charge:            '
                !> if particle charge = -1 than it is electron
                particle_charge = particle_charge * q0
                write(*, '(ES25.2)', advance='no') particle_charge
                call print_ok_mark(ierr=0)
        
            case('particle_distribution')
                write(*, '(A30, A, A)', advance='no')    '  Particle distribution:      ', &
                repeat(' ', default_length-len(trim(particle_distribution))), trim(particle_distribution)
                call print_ok_mark(ierr=0)
                if (particle_distribution == 'uniform_sphere') then
                    write(*, '(A30, ES25.2)', advance='no')    '  Sphere Center x:            ', &
                    sphere_center_x
                    if ((sphere_center_x >= x_min).and.(sphere_center_x <= x_max)) then
                        call print_ok_mark(ierr=0)
                    else
                        call print_ok_mark(ierr=1)
                    end if
                    write(*, '(A30, ES25.2)', advance='no')    '  Sphere Center y:            ', &
                    sphere_center_y
                    if ((sphere_center_y >= y_min).and.(sphere_center_y <= y_max)) then
                        call print_ok_mark(ierr=0)
                    else
                        call print_ok_mark(ierr=1)
                    end if
                    write(*, '(A30, ES25.2)', advance='no')    '  Sphere Center z:            ', &
                    sphere_center_z
                    if ((sphere_center_z >= z_min).and.(sphere_center_z <= z_max)) then
                        call print_ok_mark(ierr=0)
                    else
                        call print_ok_mark(ierr=1)
                    end if
                    write(*, '(A30, ES25.2)', advance='no')    '  Sphere Radius:              ', &
                    sphere_radius
                    if ((sphere_center_x+sphere_radius>x_max).or.(sphere_center_x-sphere_radius<x_min).or.&
                        (sphere_center_y+sphere_radius>y_max).or.(sphere_center_y-sphere_radius<y_min).or.&
                        (sphere_center_z+sphere_radius>z_max).or.(sphere_center_z-sphere_radius<z_min)) then
                        call print_ok_mark(ierr=1)
                    else
                        call print_ok_mark(ierr=0)
                    end if
                end if

            case('rotate_target')
                if (rotate_target .eqv. .true.) then
                    write(*, '(A30, L25)', advance='no')     '  Rotate target:              ', rotate_target
                    call print_ok_mark(ierr=0)
                    !> check rotate sequence
                    write(*, '(A30, A, A)', advance='no')    '  Rotate sequence:            ', &
                    repeat(' ', default_length-len(trim(rotate_sequence))), trim(rotate_sequence)
                    if ((rotate_sequence == 'xyz').or.(rotate_sequence == 'xzy').or. &
                        (rotate_sequence == 'yxz').or.(rotate_sequence == 'yzx').or. &
                        (rotate_sequence == 'zxy').or.(rotate_sequence == 'zyx').or. &
                        (rotate_sequence == 'xy').or.(rotate_sequence == 'yx').or. &
                        (rotate_sequence == 'xz').or.(rotate_sequence == 'zx').or. &
                        (rotate_sequence == 'yz').or.(rotate_sequence == 'zy').or. &
                        (rotate_sequence == 'x').or.(rotate_sequence == 'y').or. &
                        (rotate_sequence == 'z')) then
                        call print_ok_mark(ierr=0)
                    else
                        call print_ok_mark(ierr=1)
                    end if
                    !> check rotate unit
                    write(*, '(A30, A, A)', advance='no')    '  Rotate unit:                ', &
                    repeat(' ', default_length-len(trim(rotate_unit))), trim(rotate_unit)
                    if ((rotate_unit == 'degree').or.(rotate_unit == 'radian')) then
                        call print_ok_mark(ierr=0)
                    else
                        call print_ok_mark(ierr=1)
                    end if
                    !> write rotate angle
                    write(*, '(A30)', advance='no') '  Rotate about x axis:        '
                    write(*, '(F25.2)', advance='no') rotate_x
                    call print_ok_mark(ierr=0)
                    write(*, '(A30)', advance='no') '  Rotate about y axis:        '
                    write(*, '(F25.2)', advance='no') rotate_y
                    call print_ok_mark(ierr=0)
                    write(*, '(A30)', advance='no') '  Rotate about z axis:        '
                    write(*, '(F25.2)', advance='no') rotate_z
                    call print_ok_mark(ierr=0)
                end if

            case('particle_velocity_distribution')
                write(*, '(A30, A, A)', advance='no')    '  Velocity distribution:      ', &
                repeat(' ', default_length-len(trim(velocity_distribution))), trim(velocity_distribution)
                call print_ok_mark(ierr=0)

            case('temp_x')
                write(*, '(A30)', advance='no') '  Particle Temperature in x:  '
                !> temperature must >= 0
                if (particle_temp_x < 0.d0) then
                    write(*, '(ES25.2)', advance='no') particle_temp_x
                    call print_ok_mark(ierr=1)
                    stop
                else
                    write(*, '(ES25.2)', advance='no') particle_temp_x
                    call print_ok_mark(ierr=0)
                end if

            case('temp_y')
                write(*, '(A30)', advance='no') '  Particle Temperature in y:  '
                if (particle_temp_y < 0.d0) then
                    write(*, '(F25.2)', advance='no') particle_temp_y
                    call print_ok_mark(ierr=1)
                    stop
                else
                    write(*, '(F25.2)', advance='no') particle_temp_y
                    call print_ok_mark(ierr=0)
                end if

            case('temp_z')
                write(*, '(A30)', advance='no') '  Particle Temperature in z:  '
                if (particle_temp_z < 0.d0) then
                    write(*, '(F25.2)', advance='no') particle_temp_z
                    call print_ok_mark(ierr=1)
                    stop
                else
                    write(*, '(F25.2)', advance='no') particle_temp_z
                    call print_ok_mark(ierr=0)
                end if

            case('part_boundary_x')
                write(*, '(A30, A, A)', advance='no') '  Particle Boundary in x:     ', &
                repeat(' ', default_length-len(trim(part_boundary_x))), trim(part_boundary_x)
                call print_ok_mark(ierr=0)

            case('part_boundary_y')
                write(*, '(A30, A, A)', advance='no') '  Particle Boundary in y:     ', &
                repeat(' ', default_length-len(trim(part_boundary_y))), trim(part_boundary_y)
                call print_ok_mark(ierr=0)

            case('part_boundary_z')
                write(*, '(A30, A, A)', advance='no') '  Particle Boundary in z:     ', &
                repeat(' ', default_length-len(trim(part_boundary_z))), trim(part_boundary_z)
                call print_ok_mark(ierr=0)

            case('number_output_snapshots')
                write(*, '(A30)', advance='no') '  Number of output snapshots: '
                !> number of snapshot(s) must >= 0 (0 for no output)
                if (number_snapshots < 0) then
                    write(*, '(I25)', advance='no') number_snapshots
                    call print_ok_mark(ierr=1)
                    stop
                else
                    write(*, '(I25)', advance='no') number_snapshots
                    call print_ok_mark(ierr=0)
                end if

            end select
        end subroutine check_input_parameter


        subroutine print_divider(target, length)
            implicit none
            character(len=*) :: target
            integer          :: length, i, rc
            
            write(*, '(A)', advance='no') '  '
            do i = 1, length
                write(*, '(A)', advance='no') target
                rc = c_usleep(10*1000)
            end do
            write(*, '(A)') target
        end subroutine print_divider


end module helper