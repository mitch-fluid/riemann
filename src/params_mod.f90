module params_mod 
    implicit none
    private 
    public :: read_input 

    ! heat capacity ratio 
    real, public, parameter :: gamma = 1.4 
    real, public, parameter :: gm1 = gamma - 1.0 
    real, public, parameter :: gm2 = 1.0/gm1 
    real, public, parameter :: gm3 = gamma*gm2 

    ! I.C. 
    !               |
    !          2    |    1
    !      ---------|---------
    !          3    |    4
    !               |
    !
    real, public :: prim1(4), prim2(4), prim3(4), prim4(4) 

    ! dimensions and index 
    integer, public :: nghost 
    integer, public :: n, ntot 
    integer, public :: ib, ie, jb, je ! begin and end of interior pts 

    ! geometry 
    real, public, parameter :: xmin = 0.0, xmax = 1.0 
    real, public, parameter :: ymin = 0.0, ymax = 1.0 
    real, public, parameter :: xc = 0.5*(xmin + xmax)
    real, public, parameter :: yc = 0.5*(ymin + ymax) 
    real, public :: h, invh 

    ! time 
    real, public :: CFL 
    real, public :: t_final 
    integer, public :: ntmax 

contains
    subroutine read_input()
        ! local variables 
        namelist /initial/  prim1, prim2, prim3, prim4 
        namelist /spatial/  n, nghost
        namelist /temporal/ CFL, t_final, ntmax
        character(len=100) :: file_name 
        integer :: file_unit 

        ! work area 
        call get_command_argument(1, file_name) 
        open(newunit=file_unit, file=trim(file_name))
        read(file_unit, nml=initial)
        read(file_unit, nml=spatial)
        read(file_unit, nml=temporal)
        close(file_unit)

        ! print the problem info 
        call banner 

        ! index 
        ntot = n + 2*nghost 
        ib = nghost+1; ie = nghost + n 
        jb = nghost+1; je = nghost + n 

    end subroutine read_input 

    subroutine banner()
        write(*, '(a)') '!--------------------------------------------------------------------------------'
        write(*, '(a)') '!                            2D Riemann Problem                                  '
        write(*, '(a)') '!--------------------------------------------------------------------------------'
        write(*, '(a)') '!  I.C. : [rho u v p]^T                                                          '
        write(*, '(a)') '!'
        write(*, '(a)') '!        |-----------------------|-----------------------|'
        write(*, '(a, f7.4, a, f7.4,a)') '!        |        ', prim2(1), '        |        ', prim1(1), '        |'
        write(*, '(a, f7.4, a, f7.4,a)') '!        |        ', prim2(2), '        |        ', prim1(2), '        |'
        write(*, '(a, f7.4, a, f7.4,a)') '!        |        ', prim2(3), '        |        ', prim1(3), '        |'
        write(*, '(a, f7.4, a, f7.4,a)') '!        |        ', prim2(4), '        |        ', prim1(4), '        |'
        write(*, '(a)') '!        |-----------------------|-----------------------|'
        write(*, '(a, f7.4, a, f7.4,a)') '!        |        ', prim3(1), '        |        ', prim4(1), '        |'
        write(*, '(a, f7.4, a, f7.4,a)') '!        |        ', prim3(2), '        |        ', prim4(2), '        |'
        write(*, '(a, f7.4, a, f7.4,a)') '!        |        ', prim3(3), '        |        ', prim4(3), '        |'
        write(*, '(a, f7.4, a, f7.4,a)') '!        |        ', prim3(4), '        |        ', prim4(4), '        |'        
        write(*, '(a)') '!        |-----------------------|-----------------------|'
        write(*, '(a)') '!'
        write(*, '(a)') '!--------------------------------------------------------------------------------'        
        write(*, '(a, f15.7)') '!    Run until ', t_final 
        write(*, '(a)') '!--------------------------------------------------------------------------------'        
    end subroutine banner 
end module params_mod 