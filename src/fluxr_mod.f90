module fluxr_mod
    use params_mod, only: invh 
    use fhat_mod 
    implicit none
    
contains
    subroutine ffluxr(q, ql, qr, f, res, ib, ie, jb, je)
        real, intent(in) :: q(:,:,:) 
        real, intent(inout) :: f(:,:), ql(:,:), qr(:,:)
        real, intent(out) :: res(:,:,:) 
        integer, intent(in) :: ib, ie, jb, je 

        ! local variables 
        integer :: i, j 
        integer :: ib1 

        ! work area 
        ib1 = ib - 1
        do j = jb, je 
            ! copy data into the line 
            do i = ib1, ie 
                ! density 
                ql(i,1) = q(i,j,1) 
                qr(i,1) = q(i+1,j,1) 

                ! x-momentum 
                ql(i,2) = q(i,j,2)
                qr(i,2) = q(i+1,j,2) 

                ! energy 
                ql(i,3) = q(i,j,4)
                qr(i,3) = q(i+1,j,4)

                ! y-momentum 
                ql(i,4) = q(i,j,3)
                qr(i,4) = q(i+1,j,3) 
            end do 

            call fhat(ql, qr, ib1, ie, f) 

            ! compute residual 
            do i = ib, ie 
                res(i,j,1) = invh*(f(i-1,1) - f(i,1))
                res(i,j,2) = invh*(f(i-1,2) - f(i,2))
                res(i,j,3) = invh*(f(i-1,4) - f(i,4)) 
                res(i,j,4) = invh*(f(i-1,3) - f(i,3))
            end do 
        end do 
        
    end subroutine ffluxr 

    subroutine gfluxr(q, ql, qr, f, res, ib, ie, jb, je)
        real, intent(in) :: q(:,:,:) 
        real, intent(inout) :: f(:,:), ql(:,:), qr(:,:)
        real, intent(inout) :: res(:,:,:) 
        integer, intent(in) :: ib, ie, jb, je 

        ! local variables 
        integer :: i, j 
        integer :: jb1    

        ! work area 
        jb1 = jb - 1 
        do i = ib, ie 
            ! copy data into the line 
            do j = jb1, je 
                ! density 
                ql(j,1) = q(i,j,1) 
                qr(j,1) = q(i,j+1,1) 

                ! y-momentum 
                ql(j,2) = q(i,j,3)
                qr(j,2) = q(i,j+1,3) 

                ! energy 
                ql(j,3) = q(i,j,4)
                qr(j,3) = q(i,j+1,4) 

                ! x-momentum 
                ql(j,4) = q(i,j,2)
                qr(j,4) = q(i,j+1,2) 
            end do 

            call fhat(ql, qr, jb1, je, f)
            do j = jb, je 
                res(i,j,1) = res(i,j,1) + invh*(f(j-1,1) - f(j,1)) 
                res(i,j,2) = res(i,j,2) + invh*(f(j-1,4) - f(j,4)) 
                res(i,j,3) = res(i,j,3) + invh*(f(j-1,2) - f(j,2))
                res(i,j,4) = res(i,j,4) + invh*(f(j-1,3) - f(j,3))
            end do 
        end do 

    end subroutine gfluxr 
end module fluxr_mod