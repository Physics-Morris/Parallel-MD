module list
    implicit none
    double precision :: c = 299792458d0

    !> define particle type
    type, public :: particle_type
        integer          :: id
        double precision :: mass, charge
        double precision :: pos_x, pos_y, pos_z
        double precision :: vel_x, vel_y, vel_z
    contains
        procedure :: kinetic_energy => calculate_kinetic_energy
        procedure :: gamma => calculate_gamma
        procedure :: speed => calculate_speed
        procedure :: rest_energy => calculate_rest_energy

    end type particle_type

    contains

        function calculate_speed(this) result(speed)
            class(particle_type), intent(in) :: this
            double precision :: speed
             speed = dsqrt(this%vel_x**2 + this%vel_y**2 + this%vel_z**2)
        end function calculate_speed

        function calculate_gamma(this) result(gamma)
            class(particle_type), intent(in) :: this
            double precision  ::  gamma
            gamma = 1.d0 / dsqrt(1.d0 - calculate_speed(this)**2/c**2)
        end function calculate_gamma

        function calculate_kinetic_energy(this) result(KE)
            class(particle_type), intent(in) :: this
            double precision  ::  KE
            !> relativistic kinetic energy
            KE = (calculate_gamma(this) - 1.d0) * this%mass * c**2 
        end function calculate_kinetic_energy

        function calculate_rest_energy(this) result(rest_energy)
            class(particle_type), intent(in) :: this
            double precision  ::  rest_energy
            !> rest energy
            rest_energy = this%mass * c**2
        end function calculate_rest_energy

end module list


program part_list
    use list
  implicit none

  type(particle_type), dimension(:), allocatable :: local
  integer  :: i

  allocate(local(3))
  do i = 1, 3
    local(i) % id = i
    local(i) % pos_x = dble(i) * 4.d0
    local(i) % pos_y = - dble(i) * 4.d0
  end do

  do i = 1, 3
    write(*, *) local(i) % id
    write(*, *) local(i) % pos_x
    write(*, *) local(i) % pos_y
  end do
end program part_list