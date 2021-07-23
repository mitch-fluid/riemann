module misc_mod
    use data_mod 
    use params_mod 
    implicit none
    
contains 
    subroutine compute_dt(dt) 
        real, intent(out) :: dt 

        ! local variables 
        integer :: i, j 
        real :: sxmax, symax
        real :: inv, r, u, v, p, c, ke 

        ! work area 
        sxmax = 0.0
        symax = 0.0 
        do j = jb, je 
            do i = ib, ie 
                r   = q(i,j,1)
                inv = 1.0/r 
                u   = inv*q(i,j,2)
                v   = inv*q(i,j,3)
                ke  = 0.5*r*(u*u + v*v)
                p   = gm1*(q(i,j,4) - ke) 
                c   = sqrt(gamma*p*inv)

                ! estimate of the maximum wave speed 
                sxmax = max(sxmax, abs(u) + c)
                symax = max(symax, abs(v) + c) 
            end do 
        end do 

        dt = CFL*min(h/sxmax, h/symax)

    end subroutine compute_dt
end module misc_mod 