! example.f90
program main
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    implicit none

    type :: person_type
        integer           :: id
        character(len=32) :: name
        integer           :: age
    end type person_type

    integer           :: x, y
    real              :: r(2)
    type(person_type) :: alice

    ! Set initial values.
    alice = person_type(-1, 'Jane Doe', 0)
    x = 0; y = 0
    r = [ 0.0, 0.0 ]

    ! Read from file.
    call read_namelist('namelist.nml', alice, x, y, r)

    ! Output some values.
    print '(a, i0)', 'PERSON ID:   ', alice%id
    print '(2a)',    'PERSON NAME: ', alice%name
    print '(a, i0)', 'PERSON AGE:  ', alice%age
contains
    subroutine read_namelist(file_path, person, x, y, r)
        !! Reads Namelist from given file.
        character(len=*),  intent(in)  :: file_path
        type(person_type), intent(out) :: person
        integer,           intent(out) :: x, y
        real,              intent(out) :: r(2)
        integer                        :: fu, rc

        ! Namelist definition.
        namelist /EXAMPLE/ x, y, r, person

        ! Check whether file exists.
        inquire (file=file_path, iostat=rc)

        if (rc /= 0) then
            write (stderr, '(3a)') 'Error: input file "', trim(file_path), '" does not exist.'
            return
        end if

        ! Open and read Namelist file.
        open (action='read', file=file_path, iostat=rc, newunit=fu)
        read (nml=EXAMPLE, iostat=rc, unit=fu)

        if (rc /= 0) then
            write (stderr, '(a)') 'Error: invalid Namelist format.'
        end if

        close (fu)
    end subroutine read_namelist
end program main