program riemann
    use params_mod
    use data_mod 
    use init_mod
    use misc_mod
    use fluxr_mod 
    use bc_mod 
    use io_mod 
    implicit none
    real :: dt, t 
    integer :: nt 
    real :: start, finish 

    call read_input 
    call mem_alloc

    ! initialize the flow 
    call init 

    ! time loop 
    t = 0.0 
    write(*, '(a)') '        n          dt             t'
    call cpu_time(start)
    do nt = 1, ntmax 
        call compute_dt(dt) 
        if (t + dt > t_final) dt = t_final - t 

        ! enforce bc 
        call bc(q, ib, ie, jb, je)

        ! i-split 
        call ffluxr(q, ql, qr, f, res, ib, ie, jb, je)

        ! j-split 
        call gfluxr(q, ql, qr, f, res, ib, ie, jb, je) 

        ! goto 100 

        ! explicit euler 
        q = q + dt*res 

        t = t + dt 
        write(*, '(i9, 2f15.7)') nt, dt, t
        if (t == t_final) exit 
    end do 
    call cpu_time(finish) 
    write(*, '("Time elaped: ", f10.5)') finish - start

    call write_solution(q, x, y, n, n, ib, ie, jb, je) 

    call mem_free 
end program riemann 