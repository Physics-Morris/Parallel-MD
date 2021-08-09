module constants
    implicit none
    !> some constant
    double precision, PARAMETER :: pi = 3.141592653589793238462643383279503d0
    double precision, PARAMETER :: q0 = 1.602176565d-19 ! C (+/- 3.5e-27)
    double precision, PARAMETER :: m0 = 9.10938291d-31 ! kg (+/- 4e-38)
    double precision, PARAMETER :: c  = 2.99792458d8   ! m/s^2 (exact)
    double precision, PARAMETER :: c2  = 8.987551787368176d16   ! c**2
    double precision, PARAMETER :: kb = 1.3806488d-23  ! J/K (+/- 1.3e-29)
    double precision, PARAMETER :: mu0 = 4.d-7 * pi ! N/A^2 (exact)
    ! epsilon0 = 1.0_num / mu0 / c**2 ! F/m (exact)
    double precision, PARAMETER :: epsilon0 = 8.854187817620389850536563031710750d-12
    double precision, PARAMETER :: h_planck = 6.62606957d-34 ! J s (+/- 2.9e-41)
    double precision, PARAMETER :: ev = q0 ! J
    ! Derived physical parameters used in ionisation
    ! h_bar = h_planck / 2.0_num / pi
    double precision, PARAMETER :: h_bar = 1.054571725336289397963133257349698d-34

end module constants