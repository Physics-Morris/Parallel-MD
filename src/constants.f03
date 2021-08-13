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


    !> master process id
    integer, parameter          :: master_id = 0

end module constants