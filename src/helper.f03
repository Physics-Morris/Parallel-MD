module helper
    use shared_data
    use mpi
    use constants
    implicit none

    character(len=5), dimension(12) :: vt100_control = (/'[39m','[30m','[31m', &
    '[32m','[33m','[34m','[35m','[36m','[1m ','[2m ','[4m ','[0m '/)
    integer, parameter :: term_default_colour = 1
    integer, parameter :: term_black = 2
    integer, parameter :: term_red = 3
    integer, parameter :: term_green = 4
    integer, parameter :: term_yellow = 5
    integer, parameter :: term_blue = 6
    integer, parameter :: term_magnenta = 7
    integer, parameter :: term_cyan = 8
    integer, parameter :: term_bold = 9
    integer, parameter :: term_dim = 10
    integer, parameter :: term_underline = 11
    integer, parameter :: term_reset_attributes = 12
    integer, parameter :: term_max = 12

    namelist / basics_block / &
    sim_dimension, x_min, x_max, y_min, y_max, z_min, z_max

    namelist / particles_block / &
    total_particles, particle_mass, particle_charge, particle_distribution, &
    velocity_distribution

    namelist / output_block / &
    number_snapshots

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
                write(*, *) ' -----------------------------------------------'
                call set_term_color(term_default_colour)
                write(*, '(A30, I19)')     '  Simulation dimension:       ', sim_dimension
                write(*, '(A30, ES19.2)')  '  x_min of the system:        ', x_min
                write(*, '(A30, ES19.2)')  '  x_max of the system:        ', x_max
                write(*, '(A30, ES19.2)')  '  y_min of the system:        ', y_min
                write(*, '(A30, ES19.2)')  '  y_max of the system:        ', y_max
                write(*, '(A30, ES19.2)')  '  z_min of the system:        ', z_min
                write(*, '(A30, ES19.2)')  '  z_max of the system:        ', z_max
                write(*, '(A30, I19)')     '  Total number of particles:  ', total_particles
                write(*, '(A30, ES19.2)')  '  Particle mass:              ', particle_mass
                write(*, '(A30, ES19.2)')  '  Particle charge:            ', particle_charge
                write(*, '(A30, A)')       '  Particle distribution:      ', particle_distribution
                write(*, '(A30, A)')       '  Velocity distribution:      ', velocity_distribution
                write(*, '(A30, I19)')     '  Number of output snapshots: ', number_snapshots
                call set_term_color(term_default_colour)
                write(*, *) ' -----------------------------------------------'
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
                        call get_input_file('../inp/default', ierr)
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
            call sleep(1)
            call set_term_color(term_green)
            green_light = .True.
            write(*, *) 'Successful'
            call set_term_color(term_default_colour)
        end subroutine successful


        subroutine failed
            implicit none
            call sleep(1)
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
                call set_term_color(term_dim)
                write(*, '(A, ES7.1, A)', advance='no') '(', task_end_time-task_start_time, ' sec)'
                call set_term_color(term_reset_attributes)
                call set_term_color(term_bold)
            end if
        end subroutine print_task_time

end module helper
