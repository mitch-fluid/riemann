module fluxr_mod
    use roe_mod 
    implicit none
    
contains
    subroutine ffluxr(q, ql, qr, f, res, ib, ie, jb, je)
        real, intent(in) :: q(:,:,:) 
        real, intent(inout) :: f(:,:,:), ql(:,:,:), qr(:,:,:) 
        real, intent(out) :: res(:,:,:) 
        integer, intent(in) :: ib, ie, jb, je 

        ! local variables 
        integer :: i, j, l 
        integer :: ibm1, iep1  

        ! work area 
        ! copy data 
        ibm1 = ib - 1 
        iep1 = ie + 1 

        do l = 1, 4 
            do j = jb, je 
                ql(ibm1,j,l) = q(ibm1,j,l) 
                do i = ib, iep1
                    ql(i,j,l) = q(i,j,l) 
                    qr(i-1,j,l) = q(i,j,l) 
                end do 
            end do 
        end do 
        
        call roe_i(ql, qr, f, ibm1, ie, jb, je) 
        
        ! update residual 
        do l = 1, 4 
            do j = jb, je 
                do i = ib, ie 
                    res(i,j,l) = f(i-1,j,l) - f(i,j,l) 
                end do 
            end do 
        end do 

    end subroutine ffluxr 

    subroutine gfluxr(q, ql, qr, f, res, ib, ie, jb, je)
        real, intent(in) :: q(:,:,:) 
        real, intent(inout) :: f(:,:,:), ql(:,:,:), qr(:,:,:) 
        real, intent(inout) :: res(:,:,:) 
        integer, intent(in) :: ib, ie, jb, je 

        ! local variables 
        integer :: i, j, l 
        integer :: jbm1, jep1 

        jbm1 = jb - 1
        jep1 = je + 1

        ! work area 
        do l = 1, 4
            do j = jb, jep1
                do i = ib, ie 
                    ql(i,j,l) = q(i,j,l) 
                    qr(i,j-1,l) = q(i,j,l) 
                end do 
            end do 

            do i = ib, ie 
                ql(i,jbm1,l) = q(i,jbm1,l) 
            end do 
        end do 

        call roe_j(ql, qr, f, ib, ie, jbm1, je) 

        ! update residual 
        do l = 1, 4
            do j = jb, je 
                do i = ib, ie 
                    res(i,j,l) = res(i,j,l) + f(i,j-1,l) - f(i,j,l) 
                end do 
            end do 
        end do 

    end subroutine gfluxr 
end module fluxr_mod