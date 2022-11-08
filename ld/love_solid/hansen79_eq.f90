module hansen79_eq_mod
    !Ref: Penner et al. PRD 84, 103006 (2011)
    !Reformulated to a form similar to the Newtonian formalism (see e.g., Hansen, Van Horn 1979 https://ui.adsabs.harvard.edu/abs/1979ApJ...233..253H/abstract)
    !Formulation listed in Lau et al. 2019 Eqs. (25)-(30) https://journals.aps.org/prd/abstract/10.1103/PhysRevD.99.023018
    !Note: g_tt= -eNu
    !The results do not agree with Penner's formalism. Probably some mistakes in this reformulation [29-09-2022]
    implicit none
    private
    public::h79_parameters
    public::h79_eq_fluid,h79_eq_solid
    public::h79_fluid_r0_reg,h79_solid_r0_reg
    public::get_K

    integer,parameter::dp=kind(1.d0)
    integer,parameter::wp=dp

    type h79_parameters
        real(wp)::ell,p,rho,m,nu,ga,eLam,dLam,dnu,eNu
        real(wp)::mu
        real(wp)::chic ! = d mu/d P
    end type h79_parameters

    contains

    subroutine h79_eq_fluid(h79pa,t,x,f)
        !Same as p11_eq_fluid
        !x(1:2) represents {H0, H0`}
        use globalVariables,only:pi4
        type(h79_parameters),intent(in)::h79pa
        real(wp),intent(in)::t
        real(wp),dimension(2),intent(in)::x
        real(wp),dimension(2),intent(out)::f !Derivatives of {H0, H0`}
        real(wp)::l,l1
        real(wp)::r,p,rho,m,nu,ga
        real(wp)::eLam,dLam,dnu,eNu
        real(wp),dimension(2,2)::b

        l=h79pa%ell
        l1=l*(l+1.d0)
        r=t
        p=h79pa%p; rho=h79pa%rho; m=h79pa%m; nu=h79pa%nu; ga=h79pa%ga
        eLam=h79pa%eLam; dLam=h79pa%dLam; dnu=h79pa%dnu; eNu=h79pa%eNu
!        ddnu=dnu*(dLam-2.d0/r) +pi4*eLam*(rho*2+p*6-r*dnu*(rho+p))

        b=0
        b(1,2)= 1.d0
        b(2,1)= -(-l1*eLam/r**2+pi4*eLam*(rho*5+p*9+(rho+p)/(ga*p/(rho+p)) )-dnu**2)
        b(2,2)= -2.d0/r-(dnu-dLam)/2

        f(1)=dot_product(b(1,:),x)
        f(2)=dot_product(b(2,:),x)
    end subroutine h79_eq_fluid

    subroutine h79_eq_solid(h79pa,t,x,f)
        !x(1:6) represents {W,Z1,V,Z2,H0,J}
        use globalVariables,only:pi,pi4
        type(h79_parameters),intent(in)::h79pa
        real(wp),intent(in)::t
        real(wp),dimension(6),intent(in)::x
        real(wp),dimension(6),intent(out)::f !Derivatives of {W,Z1,V,Z2,H0,J}
        real(wp)::l,l1,l2
        real(wp)::r,p,rho,m,nu,ga
        real(wp)::eLam,dLam,dnu,eNu
        real(wp)::ddnu
        real(wp)::mu
        real(wp)::cs2
        real(wp)::a1,a2,a3
        real(wp),dimension(6,6)::b
        real(wp)::K,H2
        real(wp)::dpr

        l=h79pa%ell
        l1=l*(l+1)
        l2=(l+2)*(l-1)
        r=t
        p=h79pa%p; rho=h79pa%rho; m=h79pa%m; nu=h79pa%nu; ga=h79pa%ga
        eLam=h79pa%eLam; dLam=h79pa%dLam; dnu=h79pa%dnu; eNu=h79pa%eNu
        ddnu=dnu*(dLam-2/r) +pi4*eLam*(rho*2+p*6-r*dnu*(rho+p))
        mu=h79pa%mu
        cs2 = (p/(rho+p))*ga
        a1=mu; a2=cs2*(rho+p)-mu*2/3; a3=cs2*(rho+p)+mu*4/3
        dpr=-(rho+p)/2*dnu
!
        call get_K(h79pa,t,x,K)
        H2=x(5)+32*pi*a1*x(3)

        b=0
        b(1,1)=(1-2*a2/a3-r*dLam/2)/r
        b(1,2)=-r/a3
        b(1,3)=a2/a3*l1/r

        b(2,1)=(dpr*(r*ddnu/dnu-r*dLam/2-2) - 4/r*a1/a3*(a3+2*a2))/r**2
        b(2,2)=-(r*dnu/2+4*a1/a3)/r
        b(2,3)=(dpr+2/r*a1/a3*(a3+2*a2))*l1/r**2
        b(2,4)=l1*eLam/r

        b(3,1)= -eLam/r
        b(3,3)= 2/r
        b(3,4)= -r*eLam/a1

        b(4,1)= (dpr+2/r*a1/a3*(a3+2*a2))/r**2
        b(4,2)= -a2/a3/r
        b(4,3)= -(-2*a1/r+2*a1*(1+a2/a3)*l1/r)/r**2
        b(4,4)= -(r*dLam/2+r*dnu/2+3)/r
        b(4,5)= (rho+p)/2/r

        b(5,1)=(dLam+dnu)/r**2
        b(5,3)=-pi*16*dnu*a1
        b(5,6)= 1._wp

        b(6,1)=pi*32*eLam/r**2*a1/a3*(a3+2*a2) -1.5/r**2*dnu*(dLam+dnu)
        b(6,2)=-pi*8*eLam/a3*(a3+2*a2)
        b(6,3)=-pi*8/r**2*((rho+p)*eLam*l1+2*a1/a3*(a3+2*a2)*eLam*l1+4*a1*(1-eLam)-2*a1*(r*dnu)**2 )
        b(6,4)=-pi*16*eLam*(r*dnu)
        b(6,5)=(l1*eLam+2*(eLam-1)-r*(dLam/2+5*dnu/2)+(r*dnu)**2)/r**2
        b(6,6)=(r/2*(dLam-dnu)-2)/r

        f(1)=dot_product(b(1,:),x)-a2/a3*r*K-r/2*H2
        f(2)=dot_product(b(2,:),x)-(dpr+2/r*a1/a3*(a3+2*a2))*K-dpr/2*H2
        f(3)=dot_product(b(3,:),x)
        f(4)=dot_product(b(4,:),x)+a1/a3*(a3+2*a2)*K/r
        f(5)=dot_product(b(5,:),x)
        f(6)=dot_product(b(6,:),x)+((dLam+dnu)/r+pi*16*eLam*a1/a3*(a3+2*a2))*K

        f(2)=f(2)+(rho+p)/2*f(5)

    end subroutine h79_eq_solid

    subroutine get_K(h79pa,t,x,K)
        !Algebraic Eqs. (): Input {W,Zr,V,Zt,H0,J} Output: {K}
        use globalVariables,only:pi
        type(h79_parameters),intent(in)::h79pa
        real(wp),intent(in)::t
        real(wp),dimension(6),intent(in)::x
        real(wp),intent(out)::K !Perturbed metric function
        real(wp)::l,l1,l2
        real(wp)::r,p,rho,m,nu,ga
        real(wp)::eLam,dLam,dnu,eNu
        real(wp)::mu,cs2
        real(wp)::dpr
        real(wp),dimension(6)::d

        l=h79pa%ell
        l1=l*(l+1)
        l2=(l+2)*(l-1)
        r=t
        p=h79pa%p; rho=h79pa%rho; m=h79pa%m; nu=h79pa%nu; ga=h79pa%ga
        eLam=h79pa%eLam; dLam=h79pa%dLam; dnu=h79pa%dnu; eNu=h79pa%eNu
        mu=h79pa%mu
        cs2 = (p/(rho+p))*ga
        dpr=-(rho+p)/2*dnu

        d(1)=(dnu**2+dnu*dLam+pi*16*eLam*dpr*r)/(eLam*l2)
        d(2)=-pi*16*eLam*r**2/(eLam*l2)
        d(3)=0
        d(4)=-pi*16*eLam*(2+r*dnu)*r**2/(eLam*l2)
        d(5)=(l1*eLam-2+(r*dnu)**2)/(eLam*l2)
        d(6)=r**2*dnu/(eLam*l2)

        K=dot_product(d,x)
    end subroutine get_K

    function h79_fluid_r0_reg(h79pa,r)
        use globalVariables,only:pi
        type(h79_parameters),intent(in)::h79pa
        real(wp),intent(in)::r
        real(wp)::h79_fluid_r0_reg(3)
        real(wp)::l
        real(wp)::p,rho,m,nu,ga
        real(wp)::a0,a2

        l=h79pa%ell
        p=h79pa%p; rho=h79pa%rho; m=h79pa%m; nu=h79pa%nu; ga=h79pa%ga
        a0= 1._wp
        a2= -2._wp/(l*2+3._wp)*pi*(rho*5+p*9+(rho+p)/(ga*p/(rho+p)) )

        h79_fluid_r0_reg(1)= a0*r**l*(1._wp+r**2*a2)
        h79_fluid_r0_reg(2)= a0*r**(l-1._wp)*(l+r**2*(l+2._wp)*a2)
        h79_fluid_r0_reg(3)= a0*r**l*(1._wp+((l+2._wp)*a2+ pi*8/3*(p*3+rho))/(l+2._wp))   ! need second order expansion !!
    end function h79_fluid_r0_reg

    function h79_solid_r0_reg(isol,h79pa,r)
        !3 independent solutions are constructed by choosing {K0,V0,V2}
        !Follow Lau et al. 2019 https://journals.aps.org/prd/abstract/10.1103/PhysRevD.99.023018
        use globalVariables,only:pi
        type(h79_parameters),intent(in)::h79pa
        integer,intent(in)::isol
        real(wp),intent(in)::r
        real(wp)::h79_solid_r0_reg(6)
        real(wp)::l,l2
        real(wp)::p,rho,m,nu,ga
        real(wp)::mu,cs2
        real(wp)::chic
        real(wp)::H00,H02,K0,W0,W2,V0,V2,Z20,Z22
        real(wp)::Z10,Z12,J0,J2
        real(wp)::a2,a3
        real(wp)::p2,rho2
        real(wp)::mu_2nd
        real(wp)::term1,term2,term3,term4,term0

        l=h79pa%ell
        l2=(l+2._wp)*(l-1._wp)
        p=h79pa%p; rho=h79pa%rho; m=h79pa%m; nu=h79pa%nu; ga=h79pa%ga
        mu=h79pa%mu
        cs2 = (p/(rho+p))*ga
        chic=h79pa%chic

        mu_2nd=-chic*pi*2/3*(p+rho)*(p*3+rho)/2 !PRD 106, 023012 (2022) Eq. (A4)

        K0= 0; V2= 0; V0= 0
        if (isol==1) K0 = p+rho !1._wp
        if (isol==2) V2 = 1._wp/r**2
        if (isol==3) V0 = 1._wp

        W0= l*V0
        H00= K0-pi*32*mu*V0

        term1= -(1+3*cs2)/2*(rho+p)
        term2= 4._wp/9*(l*3*p*(l*p+(-l*2+6._wp)*mu)-mu_2nd*9*(l-1._wp) &
        +rho*(pi*2*(l*3*p*(2._wp-cs2)+mu*(l*(l*3-4._wp)+9)) +pi*3*l*rho*(1._wp-cs2*2) ))
        term3= (l+1._wp)/3*(mu*(l-6._wp)+l*3*(rho+p)*cs2)  ! Finn (1990) missing (l + 1)
        term0= -(l+9._wp)/3*mu -(l+3._wp)*(rho+p)*cs2
        W2=(term1*K0 + term2*V0 + term3*V2)/(-term0) ! Negative Signs account for Notation in Andersson 2011

        Z20 = -2*mu*(l-1)*V0
        Z22 = -mu*(l*V2+W2-pi*8/3*(l-2)*rho*V0)-2*mu_2nd*(l-1)*V0

        term1=pi*4*(l+3)*p-pi*4/3*(l*(2*l+3)+3)*rho-36*pi*cs2*(p+rho)
        term2=32._wp/3*pi**2*l*(p*3*(p*(3._wp+1._wp/cs2)-mu*8) &
        +rho*(p*2*(2._wp/cs2+5._wp-cs2*3)+mu*8*(l+1._wp)+rho*(1._wp/cs2+1._wp-cs2*6) ) )
        term3=-pi*8*(l+3._wp)*(rho+p)*(1._wp+cs2*3)
        term4=pi*8*l*(l+1._wp)*(rho+p)*(1._wp+cs2*3)
        term0=l*4+6._wp
        H02=(term1*K0 + term2*V0 + term3*W2 + term4*V2)/(-term0)

        Z10=-2*mu*l*(l-1)*V0
        a2=cs2*(rho+p)-2*mu/3
        a3=cs2*(rho+p)+4*mu/3
        term1=-(8*pi/3*l*rho*a3+16*pi*(a3+2*a2)*mu)-2*mu_2nd*l*(l-1)
        term2=-(a3+2*a2)/2
        term3=-(l*a3+a3+2*a2)
        term4=l*(l+1)*a2
        Z12=term1*V0+term2*H00+term3*W2+term4*V2

        J0=l*(K0-pi*8*(rho+p+4*mu)*V0)
        p2=-pi*4/3*(rho+3*p)*(rho+p)
        rho2=p2/cs2
        J2=(l+2)*H02 + pi*8/3*(pi*16*mu*(rho+3*p)*V0-3*l*(rho2+p2)*V0-3*(rho+p)*W2)

        h79_solid_r0_reg(1) = r**l*(W0+ r**2*W2)
        h79_solid_r0_reg(2) = r**(l-2)*(Z10+ r**2*Z12)
        h79_solid_r0_reg(3) = r**l*(V0+ r**2*V2)
        h79_solid_r0_reg(4) = r**(l-2)*(Z20+ r**2*Z22)
        h79_solid_r0_reg(5) = r**l*(H00+r**2*H02)
        h79_solid_r0_reg(6) = r**(l-1)*(J0+r**2*J2)
    end function h79_solid_r0_reg

end module hansen79_eq_mod

