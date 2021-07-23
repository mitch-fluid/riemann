module fhat_mod 
    use params_mod, only: gamma, gm1, gm2, gm3 
    implicit none
    private
    public :: fhat 

    real, parameter :: eps0 = 0.2 
    
contains 
    subroutine fhat(ql, qr, ib, ie, f) 
        ! arguments 
        ! q = [rho rhou rhoE rhov]^T
        real, intent(in) :: ql(:,:), qr(:,:)
        integer, intent(in) :: ib, ie 
        real, intent(out) :: f(:,:) 

        ! local variables 
        ! left and right state 
        ! density, velocity, pressure, speed of sound and 
        ! specific total enthalpy
        real :: rl, ul, vl, pl, Hl 
        real :: rr, ur, vr, pr, Hr  

        ! Roe state 
        real :: r, u, v, H, c 
        real :: wtl, wtr ! weight of left and right states

        ! Wave strength 
        real :: dr, du, dv, dp, alpha(4)

        ! Wave speed 
        real :: evsabs(4)

        ! Entropy fix 
        real :: eps 

        ! Right eigenvectors 
        real :: Rv(4,4)

        ! dissipation terms = \sum_{i=1}^4 q(\lambda)*\alpha_i*R^{(i)}
        real :: diss(4) 

        ! utility variables (avoid divide as much as possible)
        real :: inv, ke  

        ! dimension and index 
        integer :: i  
        integer :: m, n 

        ! loop over all pts 
        do i = ib, ie 
            ! left state 
            rl  = ql(i,1)
            inv = 1.0/rl 
            ul  = inv*ql(i,2)
            vl  = inv*ql(i,4)  
            ke  = 0.5*(ul*ul + vl*vl)
            pl  = gm1*(ql(i,3) - rl*ke) 
            Hl  = gm3*pl*inv + ke 

            ! right state 
            rr  = qr(i,1)
            inv = 1.0/rr 
            ur  = inv*qr(i,2)
            vr  = inv*qr(i,4)
            ke  = 0.5*(ur*ur + vr*vr) 
            pr  = gm1*(qr(i,3) - rr*ke) 
            Hr  = gm3*pr*inv + ke 

            ! Roe state 
            inv = 1.0/rr 
            wtl = sqrt(rl*inv) 
            r   = wtl*rr 
            wtl = wtl/(wtl + 1.0)
            wtr = 1.0 - wtl 
            u   = wtl*ul + wtr*ur 
            v   = wtl*vl + wtr*vr 
            H   = wtl*Hl + wtr*Hr 
            ke  = 0.5*(u*u + v*v)
            c   = sqrt(gm1*(H - ke)) 

            ! wave strength 
            dr = rr - rl 
            du = ur - ul 
            dv = vr - vl 
            dp = pr - pl 

            inv = 1.0/(c*c) 
            alpha(1) = 0.5*inv*(dp - r*c*du) 
            alpha(2) = dr - dp*inv 
            alpha(3) = 0.5*inv*(dp + r*c*du) 
            alpha(4) = r*dv 

            ! wave speed 
            evsabs(1) = abs(u - c)
            evsabs(2) = abs(u)
            evsabs(3) = abs(u + c)
            evsabs(4) = evsabs(2) 

            ! Entropy fix by Harten and Gnoffo (NASA TP-2953)
            eps = eps0*(c + abs(u) + abs(v)) 
            inv = 0.25/eps 
            if (evsabs(1) < 2.0*eps) evsabs(1) = evsabs(1)*evsabs(1)*inv + eps 
            if (evsabs(3) < 2.0*eps) evsabs(3) = evsabs(3)*evsabs(3)*inv + eps 

            ! Right eigenvectors 
            Rv(1,1) = 1.0 
            Rv(1,2) = u - c 
            Rv(1,3) = H - u*c 
            Rv(1,4) = v 

            Rv(2,1) = 1.0 
            Rv(2,2) = u 
            Rv(2,3) = 0.5*(u*u + v*v) 
            Rv(2,4) = v 

            Rv(3,1) = 1.0 
            Rv(3,2) = u + c 
            Rv(3,3) = H + u*c 
            Rv(3,4) = v 

            Rv(4,1) = 0.0 
            Rv(4,2) = 0.0 
            Rv(4,3) = v 
            Rv(4,4) = 1.0 

            ! compute dissipation 
            diss = 0.0 
            do m = 1, 4 
                do n = 1, 4 
                    diss(m) = diss(m) + evsabs(n)*alpha(n)*Rv(n,m)  
                end do 
            end do 

            ! compute flux 
            f(i,1) = 0.5*(ql(i,2) + qr(i,2) - diss(1))
            f(i,2) = 0.5*(ql(i,2)*ul + pl + qr(i,2)*ur + pr - diss(2))
            f(i,3) = 0.5*(ql(i,2)*Hl + qr(i,2)*Hr - diss(3))
            f(i,4) = 0.5*(ql(i,2)*vl + qr(i,2)*vr - diss(4)) 

        end do 
        

    end subroutine fhat 

end module fhat_mod 