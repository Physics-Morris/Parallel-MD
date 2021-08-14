module math
    use constants

    implicit none

    contains


        function uniform_distribution(vmin, vmax)
            implicit none
            double precision :: vmin, vmax, uniform_distribution
            double precision :: a, vc
            if (vmax .le. vmin) then
                write(*, *) ' error when calling uniform distribution function'
                write(*, *) ' vmax must less than vmin'
                stop
            end if
            call random_number(a)
            vc = (vmax - vmin) / 2.d0 + vmin
            a = a - 0.5d0 + vc
            a = a * 2.d0
            uniform_distribution = a * (vmax - vmin) / 2.d0
        end function uniform_distribution


        function gaussian_distribution(mu, sig)
            implicit none
            double precision :: a, b, gaussian_distribution, sig, mu
            call random_number(a)
            call random_number(b)
            gaussian_distribution = sig * dsqrt(-2.d0 * log(a)) * dcos(2.d0 * pi * b)
            gaussian_distribution = gaussian_distribution + mu
        end function gaussian_distribution


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

end module math