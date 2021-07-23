module fluxr_mod
    use fhat_mod 
    implicit none
    
contains
    subroutine ffluxr(q, ql, qr, f, res, ib, ie, jb, je)
        real, intent(in) :: q(:,:,:) 
        real, intent(inout) :: f(:,:,:), ql(:,:,:), qr(:,:,:) 
        real, intent(out) :: res(:,:,:) 
        integer, intent(in) :: ib, ie, jb, je 

        ! local variables 
        integer :: i, j 
        integer :: ibm1, iep1  

        ! work area 
        ! copy data 
        ibm1 = ib - 1 
        iep1 = ie + 1 
        do j = jb, je
            ql(ibm1,j,1) = q(ibm1,j,1) 
            ql(ibm1,j,2) = q(ibm1,j,2)
            ql(ibm1,j,3) = q(ibm1,j,4) 
            ql(ibm1,j,4) = q(ibm1,j,3)  
            do i = ib, iep1 
                ! density 
                ql(i,j,1) = q(i,j,1) 
                qr(i-1,j,1) = q(i,j,1) 

                ! x-momentum 
                ql(i,j,2) = q(i,j,2) 
                qr(i-1,j,2) = q(i,j,2) 

                ! energy 
                ql(i,j,3) = q(i,j,4)
                qr(i-1,j,3) = q(i,j,4) 

                ! y-momentum 
                ql(i,j,4) = q(i,j,3)
                qr(i-1,j,4) = q(i,j,3) 
            end do 
        end do 
        
        call fhat(ql, qr, f, ibm1, ie, jb, je) 
        
        ! update residual 
        do j = jb, je 
            do i = ib, ie 
                res(i,j,1) = f(i-1,j,1) - f(i,j,1) 
                res(i,j,2) = f(i-1,j,2) - f(i,j,2) 
                res(i,j,3) = f(i-1,j,4) - f(i,j,4) 
                res(i,j,4) = f(i-1,j,3) - f(i,j,3) 
            end do 
        end do 

    end subroutine ffluxr 

    subroutine gfluxr(q, ql, qr, f, res, ib, ie, jb, je)
        real, intent(in) :: q(:,:,:) 
        real, intent(inout) :: f(:,:,:), ql(:,:,:), qr(:,:,:) 
        real, intent(inout) :: res(:,:,:) 
        integer, intent(in) :: ib, ie, jb, je 

        ! local variables 
        integer :: i, j 
        integer :: jbm1, jep1 

        jbm1 = jb - 1
        jep1 = je + 1

        ! work area 
        do j = jb, jep1 
            do i = ib, ie 
                ! density 
                ql(i,j,1) = q(i,j,1) 
                qr(i,j-1,1) = q(i,j,1) 

                ! y-momentum 
                ql(i,j,2) = q(i,j,3) 
                qr(i,j-1,2) = q(i,j,3)

                ! energy 
                ql(i,j,3) = q(i,j,4)
                qr(i,j-1,3) = q(i,j,4) 

                ! x-momentum 
                ql(i,j,4) = q(i,j,2)
                qr(i,j-1,4) = q(i,j,2) 
            end do 
        end do 

        do i = ib, ie 
            ql(i,jbm1,1) = q(i,jbm1,1) 
            ql(i,jbm1,2) = q(i,jbm1,3) 
            ql(i,jbm1,3) = q(i,jbm1,4) 
            ql(i,jbm1,4) = q(i,jbm1,2) 
        end do 

        call fhat(ql, qr, f, ib, ie, jbm1, je) 

        ! update residual 
        do j = jb, je 
            do i = ib, ie 
                res(i,j,1) = res(i,j,1) + f(i,j-1,1) - f(i,j,1) 
                res(i,j,2) = res(i,j,2) + f(i,j-1,4) - f(i,j,4) 
                res(i,j,3) = res(i,j,3) + f(i,j-1,2) - f(i,j,2) 
                res(i,j,4) = res(i,j,4) + f(i,j-1,3) - f(i,j,3) 
            end do 
        end do 

    end subroutine gfluxr 
end module fluxr_mod