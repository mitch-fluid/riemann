module io_mod 
    use params_mod, only: gamma, gm1, gm2, gm3 
    implicit none
    private 
    public :: write_solution 

    interface write_field
        module procedure write_field_1d 
        module procedure write_field_2d 
    end interface write_field 

contains 
    subroutine write_solution(q, x, y, ni, nj, ib, ie, jb, je)
        real, intent(inout) :: q(:,:,:), x(:), y(:) 
        integer, intent(in) :: ni, nj, ib, ie, jb, je 
        ! local variables 
        integer :: i, j, l 
        real :: inv, r, u, v, p, ke 
        character(len=:), allocatable :: field_name(:)

        ! write xdmf for visualization 
        call write_xdmf(ni, nj) 

        ! transform into primitive variables (using the original memory)
        do j = jb, je 
            do i = ib, ie 
                r   = q(i,j,1)
                inv = 1.0/r 
                u   = inv*q(i,j,2)
                v   = inv*q(i,j,3) 
                ke  = 0.5*r*(u*u + v*v) 
                p   = gm1*(q(i,j,4) - ke)
                
                q(i,j,2) = u 
                q(i,j,3) = v 
                q(i,j,4) = p 
            end do 
        end do 

        ! write geomety 
        call write_field(x(ib:ie), 'x.dat') 
        call write_field(y(jb:je), 'y.dat') 

        ! write field 
        field_name = ['r.dat', 'u.dat', 'v.dat', 'p.dat']
        do l = 1, 4 
            call write_field(q(ib:ie,jb:je,l), field_name(l))
        end do 

    end subroutine write_solution  

    subroutine write_xdmf(ni, nj)
        ! arguments 
        integer, intent(in) :: ni, nj 

        ! local variables 
        integer :: file_unit, i 
        character(len=:), allocatable :: field_name(:)

        ! work area 
        field_name = ['r', 'u', 'v', 'p']
        open(newunit=file_unit, file='solution.xdmf')

        ! write header 
        write(file_unit, '(a)') '<?xml version="1.0" ?>' 
        write(file_unit, '(a)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []> '
        write(file_unit, '(a)') '<Xdmf Version="2.0">'
        write(file_unit, '(a)') '  <Domain>'        
        write(file_unit, '(a)') '    <Grid Name="Grid" GridType="Uniform">'
        write(file_unit, '(a, 2(x, i4),a)') &
        '      <Topology TopologyType="2DRectMesh" NumberOfElements="', &
        nj, ni, '"/>' 

        ! geometry 
        write(file_unit, '(a)') '      <Geometry GeometryType="VXVY">' 
        write(file_unit, '(a, x, i4, a)') &
        '        <DataItem Dimensions="', ni, &
        '" NumberType="Float" Precision="8" Format="Binary">'
        write(file_unit, '(a)') '          x.dat'
        write(file_unit, '(a)') '        </DataItem>'
        write(file_unit, '(a, x, i4, a)') &
        '        <DataItem Dimensions="', nj, &
        '" NumberType="Float" Precision="8" Format="Binary">'
        write(file_unit, '(a)') '          y.dat'
        write(file_unit, '(a)') '        </DataItem>'        
        write(file_unit, '(a)') '      </Geometry>'
        
        ! field 
        do i = 1, size(field_name)
            write(file_unit, '(a)') &
            '      <Attribute Name="'//trim(field_name(i)) &
            //'" AttributeType="Scalar" Center="Node">'
            write(file_unit, '(a, 2(x, i4), a)') &
            '        <DataItem Dimensions="', nj, ni, &
            '" NumberType="Float" Precision="8" Format="Binary">'  
            write(file_unit, '(a)') '          '//trim(field_name(i))//'.dat'
            write(file_unit, "(a)") '        </DataItem>'
            write(file_unit, "(a)") '      </Attribute>'
        end do 
        write(file_unit, "(a)") '    </Grid>' 

        ! write footer 
        write(file_unit, "(a)") '  </Domain>'
        write(file_unit, "(a)") '</Xdmf>'
        close(file_unit) 
    end subroutine write_xdmf

    subroutine write_field_1d(data, file_name)
        real, intent(in) :: data(:) 
        character(len=*), intent(in) :: file_name 

        ! local variables 
        integer :: file_unit, record_length 

        ! work area 
        inquire(iolength=record_length) data 
        open(newunit=file_unit, file=file_name, access='direct', &
             recl=record_length, form='unformatted')
        write(file_unit, rec=1) data 
        close(file_unit)
    end subroutine write_field_1d  

    subroutine write_field_2d(data, file_name)
        real, intent(in) :: data(:,:) 
        character(len=*), intent(in) :: file_name 

        ! local variables 
        integer :: file_unit, record_length 

        ! work area 
        inquire(iolength=record_length) data 
        open(newunit=file_unit, file=file_name, access='direct', &
             recl=record_length, form='unformatted')
        write(file_unit, rec=1) data 
        close(file_unit)
    end subroutine write_field_2d      
end module io_mod 