module math
    use constants
    implicit none

    contains


        function gaurand(mu, sig)
            implicit none
            double precision :: a, b, gaurand, sig, mu
            call random_number(a)
            call random_number(b)
            gaurand = sig * dsqrt(-2.d0 * log(a)) * dcos(2.d0 * pi * b) + mu
        end function


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


end module math