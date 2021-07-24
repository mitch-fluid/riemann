module data_mod 
    use params_mod
    implicit none
    
    ! geometry 
    real, allocatable, dimension(:) :: x, y 

    ! conservative variables 
    real, allocatable, dimension(:,:,:) :: q 

    ! conservative left and right state 
    real, allocatable, dimension(:,:,:) :: ql, qr 

    ! numerical flux 
    real, allocatable, dimension(:,:,:) :: f  

    ! residual 
    real, allocatable, dimension(:,:,:) :: res 

contains
    subroutine mem_alloc()
        
        allocate(x(ntot), y(ntot))
        allocate(q(ntot, ntot, 4)) 
        allocate(ql(ntot, ntot, 4), qr(ntot, ntot, 4))
        allocate(f(ntot, ntot, 4))
        allocate(res(ntot,ntot,4))

    end subroutine mem_alloc 

    subroutine mem_free()

        if (allocated(x))   deallocate(x)
        if (allocated(y))   deallocate(y) 
        if (allocated(q))   deallocate(q) 
        if (allocated(ql))  deallocate(ql)
        if (allocated(qr))  deallocate(qr)
        if (allocated(f))   deallocate(f)
        if (allocated(res)) deallocate(res)
        
    end subroutine mem_free 

end module data_mod 