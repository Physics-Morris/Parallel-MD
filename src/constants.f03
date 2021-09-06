module constants
    implicit none
    !> some constant
    double precision, parameter :: pi = 3.141592653589793238462643383279503d0
    double precision, parameter :: q0 = 1.602176565d-19 ! c (+/- 3.5e-27)
    double precision, parameter :: m0 = 9.10938291d-31 ! kg (+/- 4e-38)
    double precision, parameter :: c  = 2.99792458d8   ! m/s^2 (exact)
    double precision, parameter :: c2  = 8.987551787368176d16   ! c**2
    double precision, parameter :: kb = 1.3806488d-23  ! j/k (+/- 1.3e-29)
    double precision, parameter :: mu0 = 4.d-7 * pi ! n/a^2 (exact)
    !> epsilon0 = 1.0_num / mu0 / c**2 ! f/m (exact)
    double precision, parameter :: epsilon0 = 8.854187817620389850536563031710750d-12
    double precision, parameter :: h_planck = 6.62606957d-34 ! j s (+/- 2.9e-41)
    double precision, parameter :: ev = q0 ! j
    !> derived physical parameters used in ionisation
    !> h_bar = h_planck / 2.0_num / pi
    double precision, parameter :: h_bar = 1.054571725336289397963133257349698d-34

    !> FWHM = 2(2ln(2)) * sigma
    double precision, parameter :: fwhm2sigma=1.d0/2.d0*dsqrt(2.d0*log(2.d0))

    !> master process id
    integer, parameter          :: master_id = 0

    !> terminal related constants
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

end module constants