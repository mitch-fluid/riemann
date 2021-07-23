module init_mod
    use data_mod
    use params_mod
    implicit none
    
contains
    subroutine init()
        ! local variables 
        integer :: i, j 
        real :: r, u, v, p 
        logical :: is_left, is_right, is_lower, is_upper 

        ! work area 
        ! geometry 
        h = (xmax - xmin)/real(n-1) 
        invh = 1.0/h 
        do i = ib, ie 
            x(i) = xmin + real(i-ib)*h 
        end do 
        y(:) = x(:) 

        ! initial conditions 
        do j = jb, je 
            do i = ib, ie 
                is_left  = x(i) < xc 
                is_right = x(i) >= xc 
                is_lower = y(j) < yc 
                is_upper = y(j) >= yc 

                if (is_right .and. is_upper) then 
                    r = prim1(1)
                    u = prim1(2)
                    v = prim1(3)
                    p = prim1(4)
                else if (is_left .and. is_upper) then 
                    r = prim2(1)
                    u = prim2(2)
                    v = prim2(3)
                    p = prim2(4)
                else if (is_left .and. is_lower) then 
                    r = prim3(1)
                    u = prim3(2)
                    v = prim3(3)
                    p = prim3(4)                    
                else if (is_right .and. is_lower) then 
                    r = prim4(1)
                    u = prim4(2)
                    v = prim4(3)
                    p = prim4(4)
                end if 
                
                q(i,j,1) = r
                q(i,j,2) = r*u 
                q(i,j,3) = r*v 
                q(i,j,4) = gm2*p + 0.5*r*(u*u + v*v)                                
            end do 
        end do 

    end subroutine init 
end module init_mod 