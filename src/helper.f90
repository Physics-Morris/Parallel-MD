module helper

    contains


        subroutine get_cmd_arg
            implicit none
            character(len=*), parameter :: VERSION = '1.0'
            character(len=32)           :: arg
            integer                     :: i
        
            do i = 1, command_argument_count()
                call get_command_argument(i, arg)
        
                select case (arg)
                    case ('-v', '--version')
                        print '(2a)', 'version ', VERSION
                        stop
        
                    case ('-h', '--help')
                        call print_help()
                        stop
        
                    case default
                        print '(2a, /)', 'unrecognised command-line option: ', arg
                        call print_help()
                        stop
                end select 
            end do
        end subroutine get_cmd_arg


        subroutine print_help
            implicit none
            print '(a, /)', 'command-line options:'
            print '(a)',    '  -v, --version'
            print '(a, /)', '  -h, --help'
        end subroutine print_help


        subroutine welcome
            implicit none
            ! character(len=5), dimension(12) :: vt100_control = (/'[39m','[30m','[31m',&
            ! '[32m','[33m','[34m','[35m','[36m','[1m ','[2m ','[4m ','[0m '/)
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
            integer, intent(in) :: controlcode
            write(*,'(a)',advance='no') achar(27) // trim(vt100_control(controlcode))
        end subroutine set_term_attr


end module helper