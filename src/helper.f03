module helper
    use variables
    character(len=5), dimension(12) :: vt100_control = (/'[39m','[30m','[31m',&
    '[32m','[33m','[34m','[35m','[36m','[1m ','[2m ','[4m ','[0m '/)
    integer, parameter :: c_term_default_colour = 1
    integer, parameter :: c_term_black = 2
    integer, parameter :: c_term_red = 3
    integer, parameter :: c_term_green = 4
    integer, parameter :: c_term_yellow = 5
    integer, parameter :: c_term_blue = 6
    integer, parameter :: c_term_magenta = 7
    integer, parameter :: c_term_cyan = 8
    integer, parameter :: c_term_bold = 9
    integer, parameter :: c_term_dim = 10
    integer, parameter :: c_term_underline = 11
    integer, parameter :: c_term_reset_attributes = 12
    integer, parameter :: c_term_max = 12

    namelist /basics/ sim_dimension

    contains


        subroutine get_input_file(loc)
            implicit none
            character(len=*)            :: loc
            logical                     :: lexist

            inquire (file=loc,exist=lexist)
            if (lexist) then
                CALL set_term_attr(6)
                write(*, *)
                write(*, *) ' Input file exist !'
                write(*, *)
                CALL set_term_attr(6)
                write(*, *) ' Try to read from input file... '
                write(*, *)
                CALL set_term_attr(12)
                call sleep(2)
            else
                CALL set_term_attr(3)
                write(*, *)
                write(*, *) ' This input file ', trim(loc), ' does not exist! '
                write(*, *)
                write(*, *) ' Stopping the program... '
                write(*, *)
                CALL set_term_attr(12)
                call sleep(2)
                stop
            end if

            open(12, file=trim(loc))
            read(12, nml=basics)
            CALL set_term_attr(5)
            write(*, *) ' Reading input file successfully !'
            CALL set_term_attr(12)
            call sleep(2)
            call print_simulation_parameter
            
        end subroutine get_input_file


        subroutine print_simulation_parameter
            implicit none
            write(*, *)
            write(*, *)
            write(*, *) ' ======================================= '
            CALL set_term_attr(7)
            write(*, *) ' Simulation Dimension is', sim_dimension
            CALL set_term_attr(12)
            write(*, *) ' ======================================= '
            write(*, *)
            write(*, *)
            CALL set_term_attr(10)
            write(*, *) ' Starting program in 3 sec ... '
            CALL set_term_attr(12)
            call sleep(3)
        end subroutine print_simulation_parameter


        subroutine get_cmd_arg
            implicit none
            character(len=*), parameter :: VERSION = '1.0'
            character(len=32)           :: arg
            integer                     :: i

            if (command_argument_count() .eq. 0) then
                CALL set_term_attr(3)
                write(*, *)
                write(*, *) ' Need to sepcify input file'
                write(*, *)
                CALL set_term_attr(12)
                call print_help()
                stop
            end if

            do i = 1, command_argument_count()
                call get_command_argument(i, arg)
                select case (arg)
                    case ('-v', '--version')
                        print '(2a)', 'version ', VERSION
                        stop
        
                    case ('-h', '--help')
                        call print_help()
                        stop

                    case ('-i', '--input')
                        call get_command_argument(i+1, arg)
                        call get_input_file(arg)
                        exit

                    case ('-t', '--time')
                        call print_time
                        stop

                    case ('-d', '--default')
                        CALL set_term_attr(3)
                        write(*, *)
                        write(*, *) ' Useing default input file... '
                        CALL set_term_attr(12)
                        call get_input_file('../inp/default.input')
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
            CALL set_term_attr(5)
            write (*, *)
            call date_and_time(DATE=date, TIME=time, ZONE=zone)
            write (*, '(a, a, "-", a, "-", a, a, a, ":", a)') ' Current date and time: ', &
            date(1:4), date(5:6), date(7:8), '  ', time(1:2), time(3:4)
            write (*, *)
            CALL set_term_attr(12)
        end subroutine print_time

        subroutine print_help
            implicit none
            print '(a, /)', ' command-line options:'
            print '(a)',    '  -i, --input            , entering input file location'
            print '(a)',    '  -d, --default          , use default input file for testing'
            print '(a)',    '  -v, --version          , check out the version of the program'
            print '(a)', '  -h, --help             , see hlep for possible command line arguments'
            print '(a, /)', '  -t, --time             , print current date and time'
        end subroutine print_help


        subroutine welcome
            implicit none
            CALL set_term_attr(c_term_bold)
            CALL set_term_attr(c_term_yellow)
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*)

            CALL set_term_attr(3)
            WRITE(*,'(A)') '        ###     ###                    #########     #########       #######'
            CALL set_term_attr(1)
            WRITE(*,'(A)') '       ####   ####                    ###########   #########     ######### '
            CALL set_term_attr(2)
            CALL set_term_attr(c_term_cyan)
            CALL set_term_attr(c_term_yellow)
            WRITE(*,'(A)') '      -----------   ---    -----     ----    ----      ---       -----      '
            CALL set_term_attr(4)
            WRITE(*,'(A)') '     ### ### ###   ####   ######    ###########       ###      ######       '
            CALL set_term_attr(5)
            WRITE(*,'(A)') '    ###     ###     ##   ###       #########         ###      #######       '
            CALL set_term_attr(6)
            WRITE(*,'(A)') '   ###     ###     ##   ###       #####             ###       ######        '
            CALL set_term_attr(7)
            WRITE(*,'(A)') '  ###     ###     ##    ###      #####          ##########     ########     '
            CALL set_term_attr(8)
            WRITE(*,'(A)') ' ###     ###     ###    #####   #####          ##########      #######      '
            WRITE(*,*)

            CALL set_term_attr(c_term_green)
            WRITE(*,*)
            CALL set_term_attr(c_term_default_colour)
            CALL set_term_attr(c_term_reset_attributes)
        end subroutine welcome


        subroutine set_term_attr(controlcode)
            implicit none
            ! integer, parameter :: c_term_default_colour = 1
            ! integer, parameter :: c_term_black = 2
            ! integer, parameter :: c_term_red = 3
            ! integer, parameter :: c_term_green = 4
            ! integer, parameter :: c_term_yellow = 5
            ! integer, parameter :: c_term_blue = 6
            ! integer, parameter :: c_term_magenta = 7
            ! integer, parameter :: c_term_cyan = 8
            ! integer, parameter :: c_term_bold = 9
            ! integer, parameter :: c_term_dim = 10
            ! integer, parameter :: c_term_underline = 11
            ! integer, parameter :: c_term_reset_attributes = 12
            ! integer, parameter :: c_term_max = 12
            integer, intent(in) :: controlcode
            write(*,'(a)',advance='no') achar(27) // trim(vt100_control(controlcode))
        end subroutine set_term_attr


end module helper