module roe_mod 
    use params_mod, only: gamma, gm1, gm2, gm3 
    implicit none
    private
    public :: roe_i, roe_j 

    real, parameter :: eps = 0.2 
    real, parameter :: inveps = 1.0/eps 
    
contains 
    subroutine roe_i(ql, qr, f, ib, ie, jb, je)  
        ! arguments 
        ! q = [rho rhou rhov rhoE]^T
        real, intent(in) :: ql(:,:,:), qr(:,:,:) 
        real, intent(out) :: f(:,:,:)  
        integer, intent(in) :: ib, ie, jb, je 

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

        ! Right eigenvectors 
        real :: Ra(4,4)

        ! dissipation terms = \sum_{n=1}^4 q(\lambda_n)*\alpha_n*(R_a)_{mn}
        real :: Da(4) 

        ! utility variables (avoid divide as much as possible)
        real :: inv, ke  

        ! dimension and index 
        integer :: m, n, i, j 

        ! loop over all pts 
        do j = jb, je 
            do i = ib, ie 
                ! left state 
                rl  = ql(i,j,1)
                inv = 1.0/rl 
                ul  = inv*ql(i,j,2)
                vl  = inv*ql(i,j,3)  
                ke  = 0.5*(ul*ul + vl*vl)
                pl  = gm1*(ql(i,j,4) - rl*ke) 
                Hl  = gm3*pl*inv + ke 

                ! right state 
                rr  = qr(i,j,1)
                inv = 1.0/rr 
                ur  = inv*qr(i,j,2)
                vr  = inv*qr(i,j,3)
                ke  = 0.5*(ur*ur + vr*vr) 
                pr  = gm1*(qr(i,j,4) - rr*ke) 
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
                if (evsabs(1) < eps) evsabs(1) = 0.5*(evsabs(1)*evsabs(1)*inveps + eps)
                if (evsabs(3) < eps) evsabs(3) = 0.5*(evsabs(3)*evsabs(3)*inveps + eps) 
                
                ! Right eigenvectors 
                Ra(1,1) = 1.0 
                Ra(1,2) = u - c 
                Ra(1,3) = v 
                Ra(1,4) = H - u*c 

                Ra(2,1) = 1.0 
                Ra(2,2) = u 
                Ra(2,3) = v 
                Ra(2,4) = 0.5*(u*u + v*v) 

                Ra(3,1) = 1.0 
                Ra(3,2) = u + c 
                Ra(3,3) = v 
                Ra(3,4) = H + u*c 

                Ra(4,1) = 0.0 
                Ra(4,2) = 0.0 
                Ra(4,3) = 1.0 
                Ra(4,4) = v 

                ! compute dissipation 
                Da = 0.0 
                do m = 1, 4 
                    do n = 1, 4 
                        Da(m) = Da(m) + evsabs(n)*alpha(n)*Ra(n,m)  
                    end do 
                end do 

                ! compute flux 
                f(i,j,1) = 0.5*(ql(i,j,2) + qr(i,j,2) - Da(1))
                f(i,j,2) = 0.5*(ql(i,j,2)*ul + pl + qr(i,j,2)*ur + pr - Da(2))
                f(i,j,3) = 0.5*(ql(i,j,2)*vl + qr(i,j,2)*vr - Da(3)) 
                f(i,j,4) = 0.5*(ql(i,j,2)*Hl + qr(i,j,2)*Hr - Da(4))
                
            end do 
        end do 
        
        

    end subroutine roe_i 

    subroutine roe_j(ql, qr, f, ib, ie, jb, je)  
        ! arguments 
        ! q = [rho rhou rhov rhoE]^T
        real, intent(in) :: ql(:,:,:), qr(:,:,:) 
        real, intent(out) :: f(:,:,:)  
        integer, intent(in) :: ib, ie, jb, je 

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
        real :: dr, du, dv, dp, beta(4)

        ! Wave speed 
        real :: evsabs(4)

        ! Right eigenvectors 
        real :: Rb(4,4)

        ! dissipation terms = \sum_{n=1}^4 q(\lambda_n)*\alpha_n*(R_a)_{mn}
        real :: Db(4) 

        ! utility variables (avoid divide as much as possible)
        real :: inv, ke  

        ! dimension and index 
        integer :: m, n, i, j 

        ! loop over all pts 
        do j = jb, je 
            do i = ib, ie 
                ! left state 
                rl  = ql(i,j,1)
                inv = 1.0/rl 
                ul  = inv*ql(i,j,2)
                vl  = inv*ql(i,j,3)  
                ke  = 0.5*(ul*ul + vl*vl)
                pl  = gm1*(ql(i,j,4) - rl*ke) 
                Hl  = gm3*pl*inv + ke 

                ! right state 
                rr  = qr(i,j,1)
                inv = 1.0/rr 
                ur  = inv*qr(i,j,2)
                vr  = inv*qr(i,j,3)
                ke  = 0.5*(ur*ur + vr*vr) 
                pr  = gm1*(qr(i,j,4) - rr*ke) 
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
                beta(1) = 0.5*inv*(dp - r*c*dv) 
                beta(2) = dr - dp*inv 
                beta(3) = 0.5*inv*(dp + r*c*dv)  
                beta(4) = -r*du 

                ! wave speed 
                evsabs(1) = abs(v - c)
                evsabs(2) = abs(v)
                evsabs(3) = abs(v + c)
                evsabs(4) = evsabs(2) 

                ! Entropy fix by Harten and Gnoffo (NASA TP-2953)
                if (evsabs(1) < eps) evsabs(1) = 0.5*(evsabs(1)*evsabs(1)*inveps + eps)
                if (evsabs(3) < eps) evsabs(3) = 0.5*(evsabs(3)*evsabs(3)*inveps + eps)

                ! Right eigenvectors 
                Rb(1,1) = 1.0 
                Rb(1,2) = u 
                Rb(1,3) = v - c 
                Rb(1,4) = H - v*c 

                Rb(2,1) = 1.0 
                Rb(2,2) = u 
                Rb(2,3) = v 
                Rb(2,4) = 0.5*(u*u + v*v) 

                Rb(3,1) = 1.0 
                Rb(3,2) = u 
                Rb(3,3) = v + c 
                Rb(3,4) = H + v*c 

                Rb(4,1) = 0.0 
                Rb(4,2) = -1.0
                Rb(4,3) = 0.0 
                Rb(4,4) = -u  

                ! compute dissipation 
                Db = 0.0 
                do m = 1, 4 
                    do n = 1, 4 
                        Db(m) = Db(m) + evsabs(n)*beta(n)*Rb(n,m)  
                    end do 
                end do 

                ! compute flux 
                f(i,j,1) = 0.5*(ql(i,j,3) + qr(i,j,3) - Db(1))
                f(i,j,2) = 0.5*(ql(i,j,3)*ul + qr(i,j,3)*ur - Db(2))
                f(i,j,3) = 0.5*(ql(i,j,3)*vl + pl + qr(i,j,3)*vr + pr - Db(3)) 
                f(i,j,4) = 0.5*(ql(i,j,3)*Hl + qr(i,j,3)*Hr - Db(4)) 
                
            end do 
        end do 
        
        

    end subroutine roe_j 

end module roe_mod 