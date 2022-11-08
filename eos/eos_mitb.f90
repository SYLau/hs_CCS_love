module eos_mitb_mod
    !Mark Alford (2005) arxiv: 0411016 Eq. (1.1)
    !Modified MIT Bag Model to describe colorsuperconducting quark matter
    !P, rho in [MeV]^4 unit
    !P= fa a4 mu^4 - fa a2 mu^2 - beff
    !rho= 3 fa a4 mu^4 - fa a2 mu^2 + beff
    implicit none
    private
    public::eos_p,eos_rho,eos_ga_r
    public::eos_mu,eos_chic
    public::eos_nq
    public::mitpa

    type mit_par
        real(8)::a4,a2,beff,gap
    end type mit_par

    type(mit_par)::mitpa

    real(8),parameter::e_mks=1.60217657d-19
    real(8),parameter::hbar=1.0545718d-27
    real(8),parameter::c=2.99792458d10
    real(8),parameter::pi=3.14159265358979323846d0

    !Unit conversion:
    !rho (cgs) -> rho [(MeV)^4]
    !hbar c/L / 1E7 (cgs) -> rho (mks)
    !m c^2 / 1E7 -> E (mks)
    !E /1E6 /e (mks) -> E (MeV)

    contains

    function eos_p(rho)
        real(8),intent(in)::rho
        real(8)::eos_p
        real(8)::a4,a2,beff
        real(8)::unitconv
        real(8)::fa,muq2 !dimensionless factor; Quark chemical potential squared
        real(8)::rho_MeV4,p_MeV4

        a4=mitpa%a4; a2=mitpa%a2; beff=mitpa%beff
        unitconv=(hbar*c)**3/(1.d13*e_mks)**4*c**2 !g cm^{-3} to (MeV)^4
        fa=3.d0/4/pi**2

        rho_MeV4=rho*unitconv
        if (a2**2-12.d0*a4*(beff-rho_MeV4)/fa < 0) then
            print*,'err: eos_p Discrim < 0'
            stop
        end if
        muq2=a2/a4/6+sqrt(a2**2-12.d0*a4*(beff-rho_MeV4)/fa)/6/a4
        p_MeV4=fa*a4*muq2**2-fa*a2*muq2-beff

        eos_p =p_MeV4/unitconv*c**2
    end function eos_p

    function eos_rho(p)
        real(8),intent(in)::p
        real(8)::eos_rho
        real(8)::a4,a2,beff
        real(8)::unitconv
        real(8)::fa,muq2
        real(8)::rho_MeV4,p_MeV4

        a4=mitpa%a4; a2=mitpa%a2; beff=mitpa%beff
        unitconv=(hbar*c)**3/(1.d13*e_mks)**4 !dyn cm^{-2} to (MeV)^4
        fa=3.d0/4/pi**2

        p_MeV4=p*unitconv
        muq2=a2/a4/2+sqrt(a2**2+4.d0*a4*(beff+p_MeV4)/fa)/2/a4
        rho_MeV4=3.d0*fa*a4*muq2**2-fa*a2*muq2+beff

        eos_rho=rho_MeV4/unitconv/c**2
    end function eos_rho

    function eos_ga_r(p)
        real(8),intent(in)::p
        real(8)::eos_ga_r
        real(8)::a4,a2,beff
        real(8)::unitconv
        real(8)::fa,muq2
        real(8)::rho_MeV4,p_MeV4

        a4=mitpa%a4; a2=mitpa%a2; beff=mitpa%beff
        unitconv=(hbar*c)**3/(1.d13*e_mks)**4 !dyn cm^{-2} to (MeV)^4
        fa=3.d0/4/pi**2

        p_MeV4=p*unitconv
        muq2=a2/a4/2+sqrt(a2**2+4.d0*a4*(beff+p_MeV4)/fa)/2/a4
        rho_MeV4=3.d0*fa*a4*muq2**2-fa*a2*muq2+beff

        eos_ga_r = (rho_MeV4+p_MeV4)/p_MeV4*(a4*2*muq2-a2)/(a4*6*muq2-a2)
    end function eos_ga_r
    !gamma relativistic --- look it

    function eos_nq(p)
        real(8),intent(in)::p
        real(8)::eos_nq
        real(8)::a4,a2,beff
        real(8)::unitconv
        real(8)::fa,muq2
        real(8)::nb_MeV4,p_MeV4

        a4=mitpa%a4; a2=mitpa%a2; beff=mitpa%beff
        unitconv=(hbar*c)**3/(1.d13*e_mks)**4 !dyn cm^{-2} to (MeV)^4
        fa=3.d0/4/pi**2
!Gr= -fa a4 mu^4 + fa a2 mu^2 + beff
        p_MeV4=p*unitconv
        muq2=a2/a4/2+sqrt(a2**2+4.d0*a4*(beff+p_MeV4)/fa)/2/a4
        nb_MeV4=4*fa*a4*muq2**(3.d0/2)-2*fa*a2*sqrt(muq2)

        eos_nq=nb_MeV4/((hbar*c)**3/(1.d13*e_mks)**3)
    end function eos_nq
    !number density of quarks

    !Shear modulus of crystalline colorsuperconducting matter Mannarelli (2007)
    elemental function eos_mu(p)
        real(8),intent(in)::p
        real(8)::eos_mu
        real(8)::a4,a2,beff,gap
        real(8)::unitconv
        real(8)::fa,muq2
        real(8)::p_MeV4

        a4=mitpa%a4; a2=mitpa%a2; beff=mitpa%beff; gap=mitpa%gap
        unitconv=(hbar*c)**3/(1.d13*e_mks)**4 !dyn cm^{-2} to (MeV)^4
        fa=3.d0/4/pi**2

        p_MeV4=p*unitconv
        muq2=a2/a4/2+sqrt(a2**2+4.d0*a4*(beff+p_MeV4)/fa)/2/a4

        eos_mu=2.47d0* (gap/10)**2 * muq2/(400.d0)**2* (1.d13*e_mks*1.d39) ! convert unit from MeV/fm^3 to cgs
    end function eos_mu

    elemental function eos_chic(p)
        real(8),intent(in)::p
        real(8)::eos_chic
        real(8)::a4,a2,beff,gap
        real(8)::unitconv
        real(8)::fa,muq2
        real(8)::p_MeV4,dmu_MeV4,dp_MeV4

        a4=mitpa%a4; a2=mitpa%a2; beff=mitpa%beff; gap=mitpa%gap
        unitconv=(hbar*c)**3/(1.d13*e_mks)**4 !dyn cm^{-2} to (MeV)^4
        fa=3.d0/4/pi**2

        p_MeV4=p*unitconv
        muq2=a2/a4/2+sqrt(a2**2+4.d0*a4*(beff+p_MeV4)/fa)/2/a4
        dp_MeV4=2.d0*fa*muq2-fa*a2 !dp/d mu_q

        dmu_MeV4=2.47d0* (gap/10)**2/(400.d0)**2* (1.d13*e_mks*1.d39) *unitconv ! d mu/d mu_q; convert unit from MeV/fm^3 to (MeV)^4
        eos_chic=dmu_MeV4/dp_MeV4
    end function eos_chic

end module eos_mitb_mod

