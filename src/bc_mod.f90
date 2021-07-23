module bc_mod 
    implicit none
    
contains 
    subroutine bc(q, ib, ie, jb, je)
        real, intent(inout) :: q(:,:,:) 
        integer, intent(in) :: ib, ie, jb, je 
        ! local variables 
        integer :: i, j, l 

        ! work area 
        ! simple extrapolate 
        do l = 1, 4
            do j = jb, je 
                ! left 
                q(:ib-1,j,l) = q(ib,j,l)

                ! right 
                q(ie+1:,j,l) = q(ie,j,l) 
            end do 

            do i = ib, ie 
                ! lower 
                q(i,:jb-1,l) = q(i,jb,l)

                ! upper 
                q(i,je+1:,l) = q(i,je,l) 
            end do 
        end do 

    end subroutine bc 
end module bc_mod 