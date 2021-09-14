module math
    use constants

    implicit none

    contains


        function uniform_distribution(vmin, vmax)
            implicit none
            double precision :: vmin, vmax, uniform_distribution
            double precision :: a, width, vc
            if (vmax .le. vmin) then
                write(*, *) ' error when calling uniform distribution function'
                write(*, *) ' vmax must less than vmin'
                stop
            end if
            call random_number(a)
            width = vmax - vmin
            vc = (vmax + vmin) / 2.d0
            a = (a - 0.5d0) * width
            uniform_distribution = a + vc
        end function uniform_distribution


        !> create random gaussian distribution from Box-Muller transform
        function gaussian_distribution(mu, sig)
            implicit none
            double precision :: a, b, gaussian_distribution, sig, mu
            call random_number(a)
            call random_number(b)
            gaussian_distribution = sig * dsqrt(-2.d0 * log(a)) * dcos(2.d0 * pi * b)
            gaussian_distribution = gaussian_distribution + mu
        end function gaussian_distribution


        !> Maxwell-Boltzmann distribution has mean 0 and variance 1 with coef (kb*T/m)**0.5
        function maxwell_boltzmann(T, mass)
            implicit none
            double precision :: maxwell_boltzmann, T, mass
            maxwell_boltzmann = dsqrt(kB * T / mass) * gaussian_distribution(0.d0, 1.d0)
        end function maxwell_boltzmann


        function cross(a, b)
            implicit none
            double precision :: a(3), b(3), cross(3)
            cross(1) = a(2) * b(3) - a(3) * b(2)
            cross(2) = a(3) * b(1) - a(1) * b(3)
            cross(3) = a(1) * b(2) - a(2) * b(1)
        end function cross


        function dot(a, b)
            implicit none
            double precision :: a(3), b(3), dot
            dot = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)
        end function dot
    
    
        function mag(A)
            implicit none
            double precision :: mag, A(3) 
            mag = dsqrt(dot(A, A))
        end function
    

        function hat(A)
            implicit none
            double precision :: hat(3), A(3)
            hat = A / mag(A)
        end function


        function degree2radian(degree)
            implicit none
            double precision :: degree2radian, degree
            degree2radian = degree / 180.d0 * pi
        end function degree2radian


        function radian2degree(radian)
            implicit none
            double precision :: radian2degree, radian
            radian2degree = radian * 180.d0 / pi
        end function radian2degree


        !> rotate about x, y, z axis
        function rotate_3d(x, y, z, angle, rotate_about, unit)
            implicit none
            double precision :: rotate_3d(3), x, y, z, angle
            character(len=*) :: unit, rotate_about
            select case(rotate_about)
            !> about x-axis
            !> x'=x, y'=ycos()-zsin(), z'=ysin()+zcos()
            case ('x')
                if (unit == 'degree') then
                    rotate_3d(1) = x
                    rotate_3d(2) = y * dcos(degree2radian(angle)) - z * dsin(degree2radian(angle))
                    rotate_3d(3) = y * dsin(degree2radian(angle)) + z * dcos(degree2radian(angle))
                else if (unit == 'radian') then
                    rotate_3d(1) = x
                    rotate_3d(2) = y * dcos(angle) - z * dsin(angle)
                    rotate_3d(3) = y * dsin(angle) + z * dcos(angle)
                else
                    write(*, *) 'unrecognized option for unit'
                    stop
                end if

            !> about y-axis
            !> x'=xcos()+zsin(), y'=y, z'=-xsin()+zcos()
            case ('y')
                if (unit == 'degree') then
                    rotate_3d(1) = x * dcos(degree2radian(angle)) + z * dsin(degree2radian(angle))
                    rotate_3d(2) = y
                    rotate_3d(3) = -x * dsin(degree2radian(angle)) + z * dcos(degree2radian(angle))
                else if (unit == 'radian') then
                    rotate_3d(1) = x * dcos(angle) + z * dsin(angle)
                    rotate_3d(2) = y
                    rotate_3d(3) = -x * dsin(angle) + z * dcos(angle)
                else
                    write(*, *) 'unrecognized option for unit'
                    stop
                end if

            !> about z-axis
            !> x'=xcos()-ysin(), y'=xsin()+ycos(), z'=z
            case ('z')
                if (unit == 'degree') then
                    rotate_3d(1) = x * dcos(degree2radian(angle)) - y * dsin(degree2radian(angle))
                    rotate_3d(2) = x * dsin(degree2radian(angle)) + y * dcos(degree2radian(angle))
                    rotate_3d(3) = z
                else if (unit == 'radian') then
                    rotate_3d(1) = x * dcos(angle) - y * dsin(angle)
                    rotate_3d(2) = x * dsin(angle) + y * dcos(angle)
                    rotate_3d(3) = z
                else
                    write(*, *) 'unrecognized option for unit'
                    stop
                end if

            case default
                write(*, *) 'unrecognized option for rotate_about'
                stop

            end select

        end function rotate_3d

end module math