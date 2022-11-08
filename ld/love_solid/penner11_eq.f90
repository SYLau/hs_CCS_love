module penner11_eq_mod
    !Ref: Penner et al. PRD 84, 103006 (2011)
    !Original reference contains typo
    !Revised formalism": PRD 101, 103025 (2020)
    !Solid core regular solutions are derived from Finn MNRAS 245 82 (1990)
    !More recently derived form of the regular solutions in PRD 106, 023012 (2022) is not employed here (Haven't verified)
    !Note: g_tt= -eNu
    implicit none
    private
    public::p11_parameters
    public::p11_eq_fluid,p11_eq_solid
    public::p11_fluid_r0_reg,p11_solid_r0_reg
    public::get_euP_T1

    integer,parameter::dp=kind(1.d0)
    integer,parameter::wp=dp

    type p11_parameters
        real(wp)::ell,p,rho,m,nu,ga,eLam,dLam,dnu,eNu
        real(wp)::mu
        real(wp)::chic ! = d mu/d P
    end type p11_parameters

    contains

    subroutine p11_eq_fluid(p11pa,t,x,f)
        !x(1:2) represents {H0, H0`}
        use globalVariables,only:pi4
        type(p11_parameters),intent(in)::p11pa
        real(wp),intent(in)::t
        real(wp),dimension(2),intent(in)::x
        real(wp),dimension(2),intent(out)::f !Derivatives of {H0, H0`}
        real(wp)::l,l1
        real(wp)::r,p,rho,m,nu,ga
        real(wp)::eLam,dLam,dnu,eNu
        real(wp),dimension(2,2)::b

        l=p11pa%ell
        l1=l*(l+1._wp)
        r=t
        p=p11pa%p; rho=p11pa%rho; m=p11pa%m; nu=p11pa%nu; ga=p11pa%ga
        eLam=p11pa%eLam; dLam=p11pa%dLam; dnu=p11pa%dnu; eNu=p11pa%eNu
!        ddnu=dnu*(dLam-2._wp/r) +pi4*eLam*(rho*2+p*6-r*dnu*(rho+p))

        b=0
        b(1,2)= 1._wp
        b(2,1)= -(-l1*eLam/r**2+pi4*eLam*(rho*5+p*9+(rho+p)/(ga*p/(rho+p)) )-dnu**2)
        b(2,2)= -2._wp/r-(dnu-dLam)/2
!        b(3,1)= dnu
!        b(3,2)= 1._wp

        f(1)=dot_product(b(1,:),x)
        f(2)=dot_product(b(2,:),x)
!        f(3)=dot_product(b(3,:),x)
    end subroutine p11_eq_fluid

    subroutine p11_eq_solid(p11pa,t,x,f)
        !x(1:6) represents {H0, H0`, K, W, V, T2}
        use globalVariables,only:pi,pi4
        type(p11_parameters),intent(in)::p11pa
        real(wp),intent(in)::t
        real(wp),dimension(6),intent(in)::x
        real(wp),dimension(6),intent(out)::f !Derivatives of {H0, H0`, K, W, V, T2}
        real(wp)::l,l1,l2
        real(wp)::r,p,rho,m,nu,ga
        real(wp)::eLam,dLam,dnu,eNu
        real(wp)::mu
        real(wp)::dmu,cs2
        real(wp),dimension(6,6)::b
        real(wp)::euP,T1

        l=p11pa%ell
        l1=l*(l+1._wp)
        l2=(l+2._wp)*(l-1._wp)
        r=t
        p=p11pa%p; rho=p11pa%rho; m=p11pa%m; nu=p11pa%nu; ga=p11pa%ga
        eLam=p11pa%eLam; dLam=p11pa%dLam; dnu=p11pa%dnu; eNu=p11pa%eNu
!        ddnu=dnu*(dLam-2._wp/r) +pi4*eLam*(rho*2+p*6-r*dnu*(rho+p))
        mu=p11pa%mu
        dmu=-(rho+p)/2*dNu*p11pa%chic
        cs2 = (p/(rho+p))*ga

        call get_euP_T1(p11pa,t,x,euP,T1)

        b=0
        b(1,2)=1._wp

        b(2,1)= (l1*eLam+(eLam-1._wp)*2-r*(dLam+dnu*3)+(r*dnu)**2)/r**2
        b(2,2)= (r*(dLam-dnu)/2 -2._wp)/r
        b(2,5)= (-pi*64*mu/r**2) * (1._wp-eLam+r*(dnu+dLam/2) - 0.25d0*(r*dnu)**2) + (-pi*16*dnu)*dmu
        b(2,6)= (-pi*16*dnu/r)
!        b(2,1:6) = b(2,1:6) + (-8._wp*pi*eLam)*(3._wp+1._wp/cs2)*PP(1:6)

        b(3,1)= dnu
        b(3,2)= 1._wp
        b(3,5)= pi*16*mu/r*(r*dnu+2._wp)
        b(3,6)= -pi*16/r

        b(4,1)= -r/2
        b(4,3)= r/2
        b(4,4)= 2._wp/r-dLam/2
        b(4,5)= -(pi*16*mu*r + l1/2/r)

        b(5,4)= -eLam/r
        b(5,5)= 2._wp/r
        b(5,6)= -1._wp/mu/r

        b(6,1)= 1._wp/16/pi*(dnu+dLam)
        b(6,5)= -mu/r*eLam*l2
        b(6,6)= -((dnu-dLam)/2 + 1._wp/r)

        f(1)=dot_product(b(1,:),x)
        f(2)=dot_product(b(2,:),x) -pi*8*eLam*(3._wp+1._wp/cs2)*euP
        f(3)=dot_product(b(3,:),x)
        f(4)=dot_product(b(4,:),x) -3._wp/4/mu/r*T1
        f(5)=dot_product(b(5,:),x)
        f(6)=dot_product(b(6,:),x) -eLam*r*euP + eLam/2/r*T1

        f(2)=f(2)-pi*16*dnu*mu*f(5) !Eq. (47) contains V^\prime term
    end subroutine p11_eq_solid

    subroutine get_euP_T1(p11pa,t,x,euP,T1)
        !Algebraic Eqs. (44) (45): Input {H0,H0`,K,W,V,T2} Output: {euP, T1}
        use globalVariables,only:pi
        type(p11_parameters),intent(in)::p11pa
        real(wp),intent(in)::t
        real(wp),dimension(6),intent(in)::x
        real(wp),intent(out)::euP,T1 !Eulerian pressure perturbation, radial traction
        real(wp)::l,l1,l2
        real(wp)::r,p,rho,m,nu,ga
        real(wp)::eLam,dLam,dnu,eNu
        real(wp)::mu,cs2
        real(wp)::coef
        real(wp),dimension(2,6)::d

        l=p11pa%ell
        l1=l*(l+1._wp)
        l2=(l+2._wp)*(l-1._wp)
        r=t
        p=p11pa%p; rho=p11pa%rho; m=p11pa%m; nu=p11pa%nu; ga=p11pa%ga
        eLam=p11pa%eLam; dLam=p11pa%dLam; dnu=p11pa%dnu; eNu=p11pa%eNu
        mu=p11pa%mu
        cs2 = (p/(rho+p))*ga

        coef=pi*16*r**2*eLam- mu*4/3*r**2/(rho+p)/cs2*(-pi*16*eLam)

        d=0
        d(1,1)= l1*eLam-2._wp+(r*dnu)**2
        d(1,2)= r**2*dnu
        d(1,3)= -l2*eLam+(-pi*16*eLam)*mu*4/3 * 1.5d0*r**2
        d(1,4)= (-pi*16*eLam)*mu*4/3* (3._wp-r*dnu/2/cs2)
        d(1,5)= pi*16*mu*(r*dnu)**2+ (-pi*16*eLam)*mu*4/3 * (-1.5d0*l1)
        d(1,6)= -pi*16*(r*dnu+2._wp)
        d(1,1:6)= d(1,1:6)/coef

        d(2,3)= 1.5d0*r**2
        d(2,4)= 3._wp-r*dnu/2/cs2
        d(2,5)= -1.5d0*l1
        d(2,1:6)= d(2,1:6)+ r**2/(rho+p)/cs2*d(1,1:6)
        d(2,1:6)= d(2,1:6)*mu*4/3

        euP=dot_product(d(1,:),x)
        T1=dot_product(d(2,:),x)
    end subroutine get_euP_T1

    function p11_fluid_r0_reg(p11pa,r)
        use globalVariables,only:pi
        type(p11_parameters),intent(in)::p11pa
        real(wp),intent(in)::r
        real(wp)::p11_fluid_r0_reg(3)
        real(wp)::l
        real(wp)::p,rho,m,nu,ga
        real(wp)::a0,a2

        l=p11pa%ell
        p=p11pa%p; rho=p11pa%rho; m=p11pa%m; nu=p11pa%nu; ga=p11pa%ga
        a0= 1._wp
        a2= -2._wp/(l*2+3._wp)*pi*(rho*5+p*9+(rho+p)/(ga*p/(rho+p)) )

        p11_fluid_r0_reg(1)= a0*r**l*(1._wp+r**2*a2)
        p11_fluid_r0_reg(2)= a0*r**(l-1._wp)*(l+r**2*(l+2._wp)*a2)
        p11_fluid_r0_reg(3)= a0*r**l*(1._wp+((l+2._wp)*a2+ pi*8/3*(p*3+rho))/(l+2._wp))   ! need second order expansion !!
    end function p11_fluid_r0_reg

    function p11_solid_r0_reg(isol,p11pa,r)
        !3 independent solutions are constructed by choosing {K0,V0,V2}
        !Solid core regular solutions are derived from Finn MNRAS 245 82 (1990)
        use globalVariables,only:pi
        type(p11_parameters),intent(in)::p11pa
        integer,intent(in)::isol
        real(wp),intent(in)::r
        real(wp)::p11_solid_r0_reg(6)
        real(wp)::l,l2
        real(wp)::p,rho,m,nu,ga
        real(wp)::mu,cs2
        real(wp)::chic
        real(wp)::H00,H02,K0,K2,W0,W2,V0,V2,T20,T22
        real(wp)::mu_2nd
        real(wp)::term1,term2,term3,term4,term5,term0

        l=p11pa%ell
        l2=(l+2._wp)*(l-1._wp)
        p=p11pa%p; rho=p11pa%rho; m=p11pa%m; nu=p11pa%nu; ga=p11pa%ga
        mu=p11pa%mu
        cs2 = (p/(rho+p))*ga
        chic=p11pa%chic

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

        T20 = -mu*((l-2._wp)*V0 + W0)
        T22 = -mu*(l*V2+W2+ pi*8/3*rho*W0)-mu_2nd*((l-2._wp)*V0+W0)


        term1=pi*4*(l+3)*p-pi*4/3*(l*(2*l+3)+3)*rho-36*pi*cs2*(p+rho)
        term2=32._wp/3*pi**2*l*(p*3*(p*(3._wp+1._wp/cs2)-mu*8) &
        +rho*(p*2*(2._wp/cs2+5._wp-cs2*3)+mu*8*(l+1._wp)+rho*(1._wp/cs2+1._wp-cs2*6) ) )
        term3=-pi*8*(l+3._wp)*(rho+p)*(1._wp+cs2*3)
        term4=pi*8*l*(l+1._wp)*(rho+p)*(1._wp+cs2*3)
        term0=l*4+6._wp
        H02=(term1*K0 + term2*V0 + term3*W2 + term4*V2)/(-term0)

        term1=-pi*8/3*(rho+p*3)
        term2=pi*32/3*(pi*4*mu*((1._wp-l)*rho+p*3) -l*3*mu_2nd )
        term3=-(l+2._wp)
        term4=-pi*16*mu
        term5=-pi*16*(l+2._wp)*mu
        term0=l+2._wp
        K2=(term1*K0+term2*V0+term3*H02+term4*W2+term5*V2)/(-term0)

        p11_solid_r0_reg(1) = r**l*(H00+r**2*H02)
        p11_solid_r0_reg(2) = r**(l-1._wp)*(l*H00+r**2*(l+2._wp)*H02)
        p11_solid_r0_reg(3) = r**l*(K0+ r**2*K2)
        p11_solid_r0_reg(4) = r**l*(W0+ r**2*W2)
        p11_solid_r0_reg(5) = r**l*(V0+ r**2*V2)
        p11_solid_r0_reg(6) = r**l*(T20+ r**2*T22)
    end function p11_solid_r0_reg

    function p11_solid_r0_reg_old(isol,p11pa,r)
        !Copied from old code, divide all mu by 2, and multiply all T20, T22 by 2
        !Contains mistakes; not regular solutions
        !Left here just in case I want to do some comparisons with the old code
        !Can use it to see the divergence of the solution if the initial conditions are wrong
        use globalVariables,only:pi
        type(p11_parameters),intent(in)::p11pa
        integer,intent(in)::isol
        real(wp),intent(in)::r
        real(wp)::p11_solid_r0_reg_old(6)
        real(wp)::l,l2
        real(wp)::p,rho,m,nu,ga
        real(wp)::mu,cs2
        real(wp)::chic
        real(wp)::H00,H02,K0,K2,W0,W2,V0,V2,T20,T22
        real(wp)::mu_2nd
        real(wp)::term1,term2,term3,term4,term5

        l=p11pa%ell
        l2=(l+2._wp)*(l-1._wp)
        p=p11pa%p; rho=p11pa%rho; m=p11pa%m; nu=p11pa%nu; ga=p11pa%ga
        mu=p11pa%mu
        cs2 = (p/(rho+p))*ga
        chic=p11pa%chic

        mu_2nd=-chic*pi*2/3*(p+rho)*(p*3+rho)/2 !PRD 106, 023012 (2022) Eq. (A4)

        K0= 0; V2= 0; V0= 0
        if (isol==1) K0 = p+rho !1._wp
        if (isol==2) V2 = 1._wp/r**2
        if (isol==3) V0 = 1._wp

        W0= l*V0
        H00= K0-pi*32*mu*V0

        term1= ((ga*3+1._wp)*p + rho)/2
        term2= -(ga*p*l+ mu/3*(l-6._wp))*(l+1._wp)  ! Finn (1990) missing (l + 1) ???
        term3= (ga*p*(l+3._wp)+ mu/3*(l+9._wp))
        term4=-4._wp/9/cs2*(pi*18*mu*cs2*(-ga*6*p + p)+l*pi*9*p*(ga*p-mu*2*cs2) -cs2*9*mu_2nd*(l-1._wp) &
            +pi*rho*(l*3*(1._wp-cs2*2)*ga*p + mu*2*cs2*(l*(l*3-4._wp)-9._wp))  ) ! Some mistakes... need to update !!! [04/11/2016] See thesis !!! [29/12/2016] See thesis
        W2=(term1*(-H00) + term2*(-V2) + term4*(-V0))/term3 ! Negative Signs account for Notation in Andersson 2011

        T20 = -mu*((l-2._wp)*V0 + W0)
        T22 = -mu*(l*V2+W2+ pi*8/3*rho*W0)-mu_2nd*((l-2._wp)*V0+W0)

        term1= l*4+6._wp
        term2= pi*4/3/cs2*(-ga*p*9+cs2*((l+6._wp)*3*p -ga*p*27 -(l*(l*2+3._wp)-6._wp)*rho) )
        term3= pi**2 *32/3/cs2**2*(l*(p*3*(ga*p+cs2*3*ga*p + cs2**2*4*mu)) &
            + l*(ga*p+cs2*ga*p-cs2**2*ga*p*6+cs2**2*mu*4)*rho &
            + mu*4*cs2*(-ga*p*9+cs2*(p*12 - ga*p*27 + rho*4)) )	! corrected 3 ga*p -> -9 ga*p; cs2 9 ga*p -> cs2 -27 ga*p 27/02/2017
        term4=pi**2*64/3 * (p*3+rho)
        term5=((1._wp+cs2*3)*pi*8*ga*p)/cs2

        H02 = (term2*(-H00)+term3*(-V0)+term4*(-T20*2)-term5*(-(l+3._wp)*W2+l*(l+1._wp)*V2))/term1

        term1= pi*8*(p*3+rho)/3/(l+2._wp)
        term2= pi*64*((l-1._wp)*(l*2+1._wp)*3/2*mu_2nd+pi*4*mu*(l*rho) )/3/l2
        term3=-pi**2*32 * (p*3 + (l*2+1._wp)*rho)/3/l2
        term4= pi*32*mu/l2
        term5= pi*8*(l+1._wp)/l2

        K2 = H02+term1*H00+term2*V0+term3*T20*2+term4*(l*W2+(l*(l+1._wp)-1._wp)*V2)+term5*T22*2

        p11_solid_r0_reg_old(1) = r**l*(H00+r**2*H02)
        p11_solid_r0_reg_old(2) = r**(l-1._wp)*(l*H00+r**2*(l+2._wp)*H02)
        p11_solid_r0_reg_old(3) = r**l*(K0+ r**2*K2)
        p11_solid_r0_reg_old(4) = r**l*(W0+ r**2*W2)
        p11_solid_r0_reg_old(5) = r**l*(V0+ r**2*V2)
        p11_solid_r0_reg_old(6) = r**l*(T20+ r**2*T22)
    end function p11_solid_r0_reg_old

end module penner11_eq_mod
