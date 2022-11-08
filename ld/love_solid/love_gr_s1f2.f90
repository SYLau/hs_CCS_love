module love_gr_s1f2_mod
    !Calculate the solution of (y_out = R H0'/H0) just outside the surface !*****Check this module in detail*****
    !s1f2: layer 1 - solid; layer 2 - fluid
    !in white dwarf soldi paper, solving 6 coupled ODE, start integration and stop at density discontinuity, then use interface condition, then get initial condition for fluid layer, and then conitnue integration onto surface
    implicit none
    private
    public::get_love_gr_s1f2
    public::get_love_gr_s1f2_h79
    public::ode_eps,ode_eps_hmin
    public::ug_steps
    public::adaptive

    public::t_s,t_f,x_s,x_f

    !Integration parameters
    logical::adaptive=.true.
    real(8)::ode_eps=1.d-6
    real(8)::ode_eps_hmin=0 !Adaptive mesh method takes infinite time if gamma is not smooth; Set this value to non-zero to use adaptive
    integer::ug_steps=3000

    real(8),allocatable,dimension(:)::t_s,t_f
    real(8),allocatable,dimension(:,:)::x_s,x_f

contains
    subroutine get_love_gr_s1f2(l,savedata,yout) ! Input ell, returns R H_0^{\prime}/ H_0 evaluated at the vacuum side of the stellar surface
        use globalVariables,only:c
        use odeSolver_mod,only:rk45ad,tp,xp,rk4ug
        use bg_mod,only:tovs=>bgs,pt_i ! Treat pt_i(1) as the solid-fluid transition
        use io_mod,only:reallocate_a
!        use gnuPlot,only:gpf_plot2D
        implicit none
        real(8),intent(in)::l
        logical,intent(in)::savedata
        real(8),intent(out)::yout
        real(8)::t,x(6),tb,xv(6,3)
        real(8)::xf(2),xfi(2)
        real(8)::xsq(6)
        integer::isave,i,isol
        real(8)::h,hmin
        real(8)::emax,eps_hmin
        integer::itmax=100, iflag
        real(8)::vac_sol(2)
        real(8)::coef(3)
        integer::nold

        emax=ode_eps
        eps_hmin=ode_eps_hmin
        loop_solid: do isol=1,3
            isave=1 !save r position for faster interpolation with get_bg
            t=tovs(1)%r
            x=bc_initial(isol)

            if (size(pt_i)<1) then
                write(*,*) 'err: love_gr_s1f2_mod: pt_i must have size >= 1'
                stop
            end if

            tb=tovs(pt_i(1))%r

            if (adaptive) then
                h=(tb-t)/500
                hmin=(tb-t)*eps_hmin
                call rk45ad(ode_solid,t,x,h,tb,itmax,emax,iflag,.true.,hmin=hmin)!?
            else
                call rk4ug(ode_solid,t,x,ug_steps,tb,.true.,xsq) !?
            endif

            xv(:,isol)=x
        enddo loop_solid

        t=tovs(pt_i(1)+1)%r
        xfi=bc_interface_s1f2(pt_i(1),t,xv,coef)
        xf=xfi

        i=2
        do
            if (size(pt_i)-i+1>0) tb=tovs(pt_i(i))%r
            if (size(pt_i)-i+1==0) tb=tovs(size(tovs))%r

            if (adaptive) then
                h=(tb-t)/500
                hmin=(tb-t)*eps_hmin
                call rk45ad(ode_fluid,t,xf,h,tb,itmax,emax,iflag,.true.,hmin=hmin)
            else
                call rk4ug(ode_fluid,t,xf,ug_steps,tb,.true.,xsq(1:2))
            endif

            if (size(pt_i)-i+1==0) exit

            t=tovs(pt_i(i)+1)%r
            xf=bc_interface_f1f2(pt_i(i),t,xf)
            i=i+1
        enddo

        t=tovs(size(tovs))%r
        vac_sol=bc_surface(t,xf)
        yout=t*vac_sol(2)/vac_sol(1)

        !integrate once more with correct BCs and save solution
        if (savedata) then
            if (allocated(t_s)) deallocate(t_s)
            if (allocated(t_f)) deallocate(t_f)
            if (allocated(x_s)) deallocate(x_s)
            if (allocated(x_f)) deallocate(x_f)
            allocate(t_s(0),x_s(6,0),t_f(0),x_f(2,0))

            isave=1 !save r position for faster interpolation with get_bg
            t=tovs(1)%r
            x=coef(1)*bc_initial(1)+coef(2)*bc_initial(2)+coef(3)*bc_initial(3)

            tb=tovs(pt_i(1))%r

            if (adaptive) then
                h=(tb-t)/500
                hmin=(tb-t)*eps_hmin
                call rk45ad(ode_solid,t,x,h,tb,itmax,emax,iflag,.true.,hmin=hmin)
            else
                call rk4ug(ode_solid,t,x,ug_steps,tb,.true.,xsq)
            endif

            !save solutions
            nold=size(t_s)
            t_s=reallocate_a(t_s,nold+size(tp))
            x_s=reallocate_a(x_s,6,nold+size(tp))
            t_s(nold+1:size(t_s))=tp(1:size(tp))
            x_s(1:6,nold+1:size(t_s))=xp(1:6,:)

            t=tovs(pt_i(1)+1)%r
            xf=xfi

            i=2
            do
                if (size(pt_i)-i+1>0) tb=tovs(pt_i(i))%r
                if (size(pt_i)-i+1==0) tb=tovs(size(tovs))%r

                if (adaptive) then
                    h=(tb-t)/500
                    hmin=(tb-t)*eps_hmin
                    call rk45ad(ode_fluid,t,xf,h,tb,itmax,emax,iflag,.true.,hmin=hmin)
                else
                    call rk4ug(ode_fluid,t,xf,ug_steps,tb,.true.,xsq(1:2))
                endif
                !save solutions
                nold=size(t_f)

                t_f=reallocate_a(t_f,nold+size(tp))
                x_f=reallocate_a(x_f,2,nold+size(tp))
                t_f(nold+1:size(t_f))=tp(1:size(tp))
                x_f(1:2,nold+1:size(t_f))=xp(1:2,:)

                if (size(pt_i)-i+1==0) exit

                t=tovs(pt_i(i)+1)%r
                xf=bc_interface_f1f2(pt_i(i),t,xf)
                i=i+1
            enddo
        endif

    contains

        !Interface connecting with ODEs and BCs
        function bc_initial(isol)
            use globalVariables,only:G,c
            use penner11_eq_mod,only:p11_parameters
            use penner11_eq_mod,only:p11_solid_r0_reg
            integer,intent(in)::isol
            type(p11_parameters)::p11pa
            real(8)::P0,rho0,nu0,ga0,mu0,chi0,r0,bc_initial(6)

            P0 = G/c**4*tovs(1)%p
            rho0 = G/c**2*tovs(1)%rho
            Nu0 = tovs(1)%nu
            ga0 = tovs(1)%ga
            mu0 = G/c**4*tovs(1)%mu
            chi0 = tovs(1)%chic
            r0 = tovs(1)%r

            p11pa%ell=l
            p11pa%p=P0
            p11pa%rho=rho0
            p11pa%nu=Nu0
            p11pa%ga=ga0
            p11pa%mu=mu0
            p11pa%chic=chi0

            bc_initial=p11_solid_r0_reg(isol,p11pa,r0)
        end function

        function bc_interface_s1f2(ipos,t,xv,coef)
            use globalVariables,only:G,c,pi
            use bg_mod,only:get_tov_metric
            use penner11_eq_mod,only:p11_parameters
            use penner11_eq_mod,only:get_euP_T1
            use linsys_mod,only:LinSys_Cramer
            integer,intent(in)::ipos
            real(8),intent(in)::t,xv(6,3) !xv: index 1 - independent solution; index 2 - {H0, H0^`, K, W, V, T2}
            real(8),intent(out)::coef(3)
            real(8)::bc_interface_s1f2(2) !Output: 2 fluid pulsation quantities {H0, H0`}; Note: K is algebraically related to these two
            type(p11_parameters)::p11pa
            real(8)::gc4,gc2
            real(8)::r
            real(8)::p,rho,m,nu,ga
            real(8)::mu
            real(8)::eLam,dLam,dnu,eNu
            real(8)::rho_fluid
            real(8),dimension(3)::euP,T1 !Eulerian pressure perturbation, radial traction
            real(8),dimension(2,2)::U
            real(8),dimension(2)::vec
            real(8)::H0,H0prime,K
            integer::i

            gc2=G/c**2
            gc4=G/c**4

            !solid side variables at the interface
            r=t
            p=tovs(ipos)%p
            rho=tovs(ipos)%rho
            m=tovs(ipos)%m
            nu=tovs(ipos)%nu
            ga=tovs(ipos)%ga
            mu=tovs(ipos)%mu
            call get_tov_metric(r,p,rho,m,nu,eLam,dLam,dnu,eNu)
            p=gc4*p; rho=gc2*rho; m=gc2*m !Unit conversion
            mu=gc4*mu

            p11pa%ell=l
            p11pa%p=p; p11pa%rho=rho; p11pa%m=m; p11pa%nu=nu; p11pa%ga=ga
            p11pa%eLam=eLam; p11pa%dLam=dLam; p11pa%dnu=dnu; p11pa%eNu=eNu
            p11pa%mu=mu
            do i=1,3
                call get_euP_T1(p11pa,r,xv(:,i),euP(i),T1(i))
            end do

            !fluid side density
            rho_fluid=gc2*tovs(ipos+1)%rho

            coef(1)=1.d0
            !Condition 1: Eq.(42), (50) of Penner 2011: r^2 (euP + dP/dr W/r) + T1 = r^2 (euP + dP/dr W/r)_{fluid}
            ! = r^2 [(rho+P)/2 *H0 + dP/dr W/r]_{fluid}
            !xv(1,:) = H0; xv(4,:) = W
            U(1,:)=r**2*euP(2:3)+T1(2:3)-r*(rho-rho_fluid)/2*dNu*xv(4,2:3)-r**2*(rho_fluid+p)/2*xv(1,2:3)
            vec(1)=-(r**2*euP(1)+T1(1)-r*(rho-rho_fluid)/2*dNu*xv(4,1)-r**2*(rho_fluid+p)/2*xv(1,1))*coef(1)

            !Condition 2: T2 = 0
            !xv(6,:) = T2
            U(2,:)=xv(6,2:3)
            vec(2)=-xv(6,1)*coef(1)

            call LinSys_Cramer(U,vec,coef(2:3))

            !Algebraic equation
            !(l(l+1)-2)eLam K =  (r^2 dNu)H0` + [l(l+1)eLam - 2 + (r dNu)^2 - r(dNu+dLam)]H0
            !dNu+dLam = 8pi r eLam (rho+P)
            ![Quick check: Fluid limit of Penner (2011) Eq. (38)]
            H0=dot_product(coef,xv(1,:))
            K=dot_product(coef,xv(3,:))
            H0prime=(eLam*(l*(l+1)-2)*K - (eLam*l*(l+1)-2+(r*dNu)**2-pi*8*eLam*r**2*(rho_fluid+p))*H0)/(r**2*dNu)

            bc_interface_s1f2(1)=H0
            bc_interface_s1f2(2)=H0prime
        end function bc_interface_s1f2

        function bc_interface_f1f2(ipos,t,x)
            use globalVariables,only:pi4,c
            integer,intent(in)::ipos
            real(8),intent(in)::t,x(2)
            real(8)::bc_interface_f1f2(2)
            real(8)::r,m_i,p_i,drho_i

            r=t
            m_i=tovs(ipos+1)%m
            p_i=tovs(ipos+1)%p
            drho_i=tovs(ipos+1)%rho-tovs(ipos)%rho
            bc_interface_f1f2(1)=x(1)
            bc_interface_f1f2(2)=x(2)+ pi4*r**2/(m_i+pi4*r**3*p_i/c**2)*drho_i*x(1) !Junction condition when there is density discontinuity
        end function bc_interface_f1f2

        function bc_surface(t,x)
            use globalVariables,only:pi4
            real(8),intent(in)::t,x(2)
            real(8)::bc_surface(2)
            real(8)::r

            r=t
            bc_surface(1)=x(1)
            bc_surface(2)=x(2)-pi4*r**2/(tovs(size(tovs))%m)*tovs(size(tovs))%rho*x(1)
        end function

        subroutine ode_fluid(t,x,f)
            use globalVariables,only:G,c,pi4
            use bg_mod,only:get_tov_metric,get_bg_r
            use penner11_eq_mod,only:p11_parameters,p11_eq_fluid
            real(8), intent(in):: t, x(:)
            real(8), intent(out):: f(:)
            real(8)::gc2,gc4,p,rho,m,nu,ga,A,eLam,dLam,dnu,eNu
            type(p11_parameters)::p11pa

            gc2=G/c**2
            gc4=G/c**4

            call get_bg_r(t,p,rho,m,nu,ga,A,isave)
            call get_tov_metric(t,p,rho,m,nu,eLam,dLam,dnu,eNu)
            p=gc4*p; rho=gc2*rho; m=gc2*m !unit conversion

            p11pa%ell=l
            p11pa%p=p; p11pa%rho=rho; p11pa%m=m; p11pa%nu=nu; p11pa%ga=ga
            p11pa%eLam=eLam; p11pa%dLam=dLam; p11pa%dnu=dnu; p11pa%eNu=eNu

            call p11_eq_fluid(p11pa,t,x,f)
        end subroutine


        subroutine ode_solid(t,x,f)
            use globalVariables,only:G,c,pi4
            use bg_mod,only:get_tov_metric,get_bg_solid_r
            use penner11_eq_mod,only:p11_parameters,p11_eq_solid
            real(8), intent(in):: t, x(:)
            real(8), intent(out):: f(:)
            real(8)::gc2,gc4,p,rho,m,nu,ga,A,mu,chic,eLam,dLam,dnu,eNu
            type(p11_parameters)::p11pa

            gc2=G/c**2
            gc4=G/c**4

            call get_bg_solid_r(t,p,rho,m,nu,ga,A,mu,chic,isave)
            call get_tov_metric(t,p,rho,m,nu,eLam,dLam,dnu,eNu)
            p=gc4*p; rho=gc2*rho; m=gc2*m; mu=gc4*mu !unit conversion

            p11pa%ell=l
            p11pa%p=p; p11pa%rho=rho; p11pa%m=m; p11pa%nu=nu; p11pa%ga=ga
            p11pa%eLam=eLam; p11pa%dLam=dLam; p11pa%dnu=dnu; p11pa%eNu=eNu
            p11pa%mu=mu; p11pa%chic=chic

            call p11_eq_solid(p11pa,t,x,f)
        end subroutine

!#######################################################################################################

    end subroutine get_love_gr_s1f2

    !alternative formalism in-progress,
    subroutine get_love_gr_s1f2_h79(l,savedata,yout) ! Input ell, returns R H_0^{\prime}/ H_0 evaluated at the vacuum side of the stellar surface
        use globalVariables,only:c
        use odeSolver_mod,only:rk45ad,tp,xp,rk4ug
        use bg_mod,only:tovs=>bgs,pt_i ! Treat pt_i(1) as the solid-fluid transition
        use io_mod,only:reallocate_a
        implicit none
        real(8),intent(in)::l
        logical,intent(in)::savedata
        real(8),intent(out)::yout
        real(8)::t,x(6),tb,xv(6,3)
        real(8)::xf(2),xfi(2)
        real(8)::xsq(6)
        integer::isave,i,isol
        real(8)::h,hmin
        real(8)::emax,eps_hmin
        integer::itmax=100, iflag
        real(8)::vac_sol(2)
        real(8)::coef(3)
        integer::nold

        emax=ode_eps
        eps_hmin=ode_eps_hmin
        loop_solid: do isol=1,3
            isave=1 !save r position for faster interpolation with get_bg
            t=tovs(1)%r
            x=bc_initial(isol)

            if (size(pt_i)<1) then
                write(*,*) 'err: love_gr_s1f2_mod: pt_i must have size >= 1'
                stop
            end if

            tb=tovs(pt_i(1))%r

            if (adaptive) then
                h=(tb-t)/500
                hmin=(tb-t)*eps_hmin
                call rk45ad(ode_solid,t,x,h,tb,itmax,emax,iflag,.true.,hmin=hmin)
            else
                call rk4ug(ode_solid,t,x,ug_steps,tb,.true.,xsq)
            endif

            xv(:,isol)=x
        enddo loop_solid

        t=tovs(pt_i(1)+1)%r
        xfi=bc_interface_s1f2(pt_i(1),t,xv,coef)
        xf=xfi

        i=2
        do
            if (size(pt_i)-i+1>0) tb=tovs(pt_i(i))%r
            if (size(pt_i)-i+1==0) tb=tovs(size(tovs))%r

            if (adaptive) then
                h=(tb-t)/500
                hmin=(tb-t)*eps_hmin
                call rk45ad(ode_fluid,t,xf,h,tb,itmax,emax,iflag,.true.,hmin=hmin)
            else
                call rk4ug(ode_fluid,t,xf,ug_steps,tb,.true.,xsq(1:2))
            endif

            if (size(pt_i)-i+1==0) exit

            t=tovs(pt_i(i)+1)%r
            xf=bc_interface_f1f2(pt_i(i),t,xf)
            i=i+1
        enddo

        t=tovs(size(tovs))%r
        vac_sol=bc_surface(t,xf)
        yout=t*vac_sol(2)/vac_sol(1)

        !integrate once more with correct BCs and save solution
        if (savedata) then
            if (allocated(t_s)) deallocate(t_s)
            if (allocated(t_f)) deallocate(t_f)
            if (allocated(x_s)) deallocate(x_s)
            if (allocated(x_f)) deallocate(x_f)
            allocate(t_s(0),x_s(6,0),t_f(0),x_f(2,0))

            isave=1 !save r position for faster interpolation with get_bg
            t=tovs(1)%r
            x=coef(1)*bc_initial(1)+coef(2)*bc_initial(2)+coef(3)*bc_initial(3)

            tb=tovs(pt_i(1))%r

            if (adaptive) then
                h=(tb-t)/500
                hmin=(tb-t)*eps_hmin
                call rk45ad(ode_solid,t,x,h,tb,itmax,emax,iflag,.true.,hmin=hmin)
            else
                call rk4ug(ode_solid,t,x,ug_steps,tb,.true.,xsq)
            endif

            !save solutions
            nold=size(t_s)
            t_s=reallocate_a(t_s,nold+size(tp))
            x_s=reallocate_a(x_s,6,nold+size(tp))
            t_s(nold+1:size(t_s))=tp(1:size(tp))
            x_s(1:6,nold+1:size(t_s))=xp(1:6,:)

            t=tovs(pt_i(1)+1)%r
            xf=xfi

            i=2
            do
                if (size(pt_i)-i+1>0) tb=tovs(pt_i(i))%r
                if (size(pt_i)-i+1==0) tb=tovs(size(tovs))%r

                if (adaptive) then
                    h=(tb-t)/500
                    hmin=(tb-t)*eps_hmin
                    call rk45ad(ode_fluid,t,xf,h,tb,itmax,emax,iflag,.true.,hmin=hmin)
                else
                    call rk4ug(ode_fluid,t,xf,ug_steps,tb,.true.,xsq(1:2))
                endif
                !save solutions
                nold=size(t_f)

                t_f=reallocate_a(t_f,nold+size(tp))
                x_f=reallocate_a(x_f,2,nold+size(tp))
                t_f(nold+1:size(t_f))=tp(1:size(tp))
                x_f(1:2,nold+1:size(t_f))=xp(1:2,:)

                if (size(pt_i)-i+1==0) exit

                t=tovs(pt_i(i)+1)%r
                xf=bc_interface_f1f2(pt_i(i),t,xf)
                i=i+1
            enddo
        endif

    contains

        !Interface connecting with ODEs and BCs
        function bc_initial(isol)
            use globalVariables,only:G,c
            use hansen79_eq_mod,only:h79_parameters
            use hansen79_eq_mod,only:h79_solid_r0_reg
            integer,intent(in)::isol
            type(h79_parameters)::h79pa
            real(8)::P0,rho0,nu0,ga0,mu0,chi0,r0,bc_initial(6)

            P0 = G/c**4*tovs(1)%p
            rho0 = G/c**2*tovs(1)%rho
            Nu0 = tovs(1)%nu
            ga0 = tovs(1)%ga
            mu0 = G/c**4*tovs(1)%mu
            chi0 = tovs(1)%chic
            r0 = tovs(1)%r

            h79pa%ell=l
            h79pa%p=P0
            h79pa%rho=rho0
            h79pa%nu=Nu0
            h79pa%ga=ga0
            h79pa%mu=mu0
            h79pa%chic=chi0

            bc_initial=h79_solid_r0_reg(isol,h79pa,r0)
        end function

        function bc_interface_s1f2(ipos,t,xv,coef)
            use globalVariables,only:G,c,pi
            use bg_mod,only:get_tov_metric
            use hansen79_eq_mod,only:h79_parameters
            use linsys_mod,only:LinSys_Cramer
            integer,intent(in)::ipos
            real(8),intent(in)::t,xv(6,3) !xv: index 1 - independent solution; index 2 - {W,Zr,V,Zt,H0,J}
            real(8),intent(out)::coef(3)
            real(8)::bc_interface_s1f2(2) !Output: 2 fluid pulsation quantities {H0, H0`}; Note: K is algebraically related to these two
            type(h79_parameters)::h79pa
            real(8)::gc4,gc2
            real(8)::r
            real(8)::p,rho,m,nu,ga
            real(8)::mu
            real(8)::eLam,dLam,dnu,eNu
            real(8)::rho_fluid
            real(8),dimension(2,2)::U
            real(8),dimension(2)::vec
            real(8)::H0,H0prime

            gc2=G/c**2
            gc4=G/c**4

            !solid side variables at the interface
            r=t
            p=tovs(ipos)%p
            rho=tovs(ipos)%rho
            m=tovs(ipos)%m
            nu=tovs(ipos)%nu
            ga=tovs(ipos)%ga
            mu=tovs(ipos)%mu
            call get_tov_metric(r,p,rho,m,nu,eLam,dLam,dnu,eNu)
            p=gc4*p; rho=gc2*rho; m=gc2*m !Unit conversion
            mu=gc4*mu

            h79pa%ell=l
            h79pa%p=p; h79pa%rho=rho; h79pa%m=m; h79pa%nu=nu; h79pa%ga=ga
            h79pa%eLam=eLam; h79pa%dLam=dLam; h79pa%dnu=dnu; h79pa%eNu=eNu
            h79pa%mu=mu

            !fluid side density
            rho_fluid=gc2*tovs(ipos+1)%rho

            coef(1)=1.d0
            !Condition 1: Z1=(rho_fluid+ p)/2 * (H0 - dnu*W/R)
            !xv(1,:) = W; xv(2,:) = Z1; xv(5,:) = H0
            U(1,:)=(xv(2,2:3)-(rho_fluid+p)/2*(xv(5,2:3)-dnu*xv(1,2:3)/r))
            vec(1)=-(xv(2,1)-(rho_fluid+p)/2*(xv(5,1)-dnu*xv(1,1)/r))*coef(1)

            !Condition 2: Z22 = 0
            !xv(4,:) = Z22
            U(2,:)=xv(4,2:3)
            vec(2)=-xv(4,1)*coef(1)

            call LinSys_Cramer(U,vec,coef(2:3))

            !Algebraic equation
            !(l(l+1)-2)eLam K =  (r^2 dNu)H0` + [l(l+1)eLam - 2 + (r dNu)^2 - r(dNu+dLam)]H0
            !dNu+dLam = 8pi r eLam (rho+P)
            ![Quick check: Fluid limit of Penner (2011) Eq. (38)]
            H0=dot_product(coef,xv(5,:))
            H0prime=dot_product(coef,xv(6,:))

            bc_interface_s1f2(1)=H0
            bc_interface_s1f2(2)=H0prime
        end function bc_interface_s1f2

        function bc_interface_f1f2(ipos,t,x)
            use globalVariables,only:pi4,c
            integer,intent(in)::ipos
            real(8),intent(in)::t,x(2)
            real(8)::bc_interface_f1f2(2)
            real(8)::r,m_i,p_i,drho_i

            r=t
            m_i=tovs(ipos+1)%m
            p_i=tovs(ipos+1)%p
            drho_i=tovs(ipos+1)%rho-tovs(ipos)%rho
            bc_interface_f1f2(1)=x(1)
            bc_interface_f1f2(2)=x(2)+ pi4*r**2/(m_i+pi4*r**3*p_i/c**2)*drho_i*x(1) !Junction condition when there is density discontinuity
        end function bc_interface_f1f2

        function bc_surface(t,x)
            use globalVariables,only:pi4
            real(8),intent(in)::t,x(2)
            real(8)::bc_surface(2)
            real(8)::r

            r=t
            bc_surface(1)=x(1)
            bc_surface(2)=x(2)-pi4*r**2/(tovs(size(tovs))%m)*tovs(size(tovs))%rho*x(1)
        end function

        subroutine ode_fluid(t,x,f)
            use globalVariables,only:G,c,pi4
            use bg_mod,only:get_tov_metric,get_bg_r
            use hansen79_eq_mod,only:h79_parameters,h79_eq_fluid
            real(8), intent(in):: t, x(:)
            real(8), intent(out):: f(:)
            real(8)::gc2,gc4,p,rho,m,nu,ga,A,eLam,dLam,dnu,eNu
            type(h79_parameters)::h79pa

            gc2=G/c**2
            gc4=G/c**4

            call get_bg_r(t,p,rho,m,nu,ga,A,isave)
            call get_tov_metric(t,p,rho,m,nu,eLam,dLam,dnu,eNu)
            p=gc4*p; rho=gc2*rho; m=gc2*m !unit conversion

            h79pa%ell=l
            h79pa%p=p; h79pa%rho=rho; h79pa%m=m; h79pa%nu=nu; h79pa%ga=ga
            h79pa%eLam=eLam; h79pa%dLam=dLam; h79pa%dnu=dnu; h79pa%eNu=eNu

            call h79_eq_fluid(h79pa,t,x,f)
        end subroutine


        subroutine ode_solid(t,x,f)
            use globalVariables,only:G,c,pi4
            use bg_mod,only:get_tov_metric,get_bg_solid_r
            use hansen79_eq_mod,only:h79_parameters,h79_eq_solid
            real(8), intent(in):: t, x(:)
            real(8), intent(out):: f(:)
            real(8)::gc2,gc4,p,rho,m,nu,ga,A,mu,chic,eLam,dLam,dnu,eNu
            type(h79_parameters)::h79pa

            gc2=G/c**2
            gc4=G/c**4

            call get_bg_solid_r(t,p,rho,m,nu,ga,A,mu,chic,isave)
            call get_tov_metric(t,p,rho,m,nu,eLam,dLam,dnu,eNu)
            p=gc4*p; rho=gc2*rho; m=gc2*m; mu=gc4*mu !unit conversion

            h79pa%ell=l
            h79pa%p=p; h79pa%rho=rho; h79pa%m=m; h79pa%nu=nu; h79pa%ga=ga
            h79pa%eLam=eLam; h79pa%dLam=dLam; h79pa%dnu=dnu; h79pa%eNu=eNu
            h79pa%mu=mu; h79pa%chic=chic

            call h79_eq_solid(h79pa,t,x,f)
        end subroutine

    end subroutine get_love_gr_s1f2_h79

end module

