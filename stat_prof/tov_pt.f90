module tov_pt_mod
    !This module solves the TOV equations when provided the EOS in the format specified by the interface below [see the subroutine tov_pt_driver]
    implicit none
    private
    public::tov_pt_bg
    public::tov_pt_driver
    public::tov_pt_get_mr,tov_pt_get_rhoc

    type tov_pt_bg
        real(8)::r,p,rho,m,nu,ga,A !Schwarzchild discriminant "A" default to be zero (barotrope)
    end type
    real(8),allocatable,private:: t_bg(:),x_bg(:,:)
    !type(tov_pt_bg),allocatable,protected,public:: tovs(:)
    !integer,allocatable,protected,public:: pt_i(:)

    real(8),parameter,private::eps_r0=1.d-7 !7.15d-7 !initial radius

contains

    subroutine tov_pt_driver(rhoc,pt,eps_Pm,eos_p,eos_rho,eos_ga_r,tovs,pt_i) !solve TOV eq for models with multiple 1st order phase transitions; store results in tovs (derived data type)
        use globalVariables,only:G,pi,c
        use nr_mod,only:polint
        implicit none
        real(8),intent(in)::rhoc,pt(:),eps_Pm
        type(tov_pt_bg),allocatable,intent(out)::tovs(:)
        integer,allocatable,intent(out):: pt_i(:)
        real(8)::r0,P0,m0,nu0,Pc,rc,R2,Pm,h0,nuR
        real(8)::err
        integer::i,j,nold
        interface
            function eos_rho(p) !input p return rho
                real(8),intent(in)::p
                real(8)::eos_rho
            end function eos_rho
            function eos_p(rho)
                real(8),intent(in)::rho
                real(8)::eos_p
            end function eos_p
            function eos_ga_r(p)
                real(8),intent(in)::p
                real(8)::eos_ga_r
            end function eos_ga_r
        end interface

        Pc= eos_p(rhoc)
        R2=sqrt(3.d0*Pc/(2.d0*pi*G*(rhoc+Pc/c**2)*(rhoc+3.d0*Pc/c**2))) !Taylor's expansion of dP/dr; see e.g. Kruger et al. 2015
        rc= eps_r0*R2

        !if (allocated(tovs)) deallocate(tovs)
        !if (allocated(pt_i)) deallocate(pt_i)
        allocate(tovs(0),pt_i(size(pt)))
        i=1
        r0=rc
        P0=Pc*(1.d0-(rc/R2)**2)
        m0=4.d0*pi/3.d0*rhoc*rc**3
        nu0=2.d0*pi/3.d0*G/c**2*(rhoc+3.d0*Pc/c**2)*rc**2
        h0=rc

        do
            if (size(pt)-i+1>0) Pm=pt(i)
            if (size(pt)-i+1==0) Pm=eps_Pm*Pc
            call integrate_tov(P0,m0,nu0,r0,Pm,h0,eos_rho)

            nold=size(tovs)
            tovs=reallocate_bg_a(tovs,size(tovs)+size(t_bg))
            if (size(pt)-i+1>0) pt_i(i)=size(tovs)
            do j=1,size(t_bg)-1
                tovs(j+nold)=tov_pt_bg(t_bg(j),x_bg(1,j),x_bg(2,j),x_bg(3,j),x_bg(4,j),eos_ga_r(x_bg(1,j)),0.d0)
            end do
            !interpolation for last data point
            tovs(size(tovs))%P=Pm
            tovs(size(tovs))%rho=eos_rho(Pm*(1.d0+eps_Pm))
            tovs(size(tovs))%ga=eos_ga_r(Pm*(1.d0+eps_Pm))

            call polint(x_bg(1,size(t_bg)-1:size(t_bg)),t_bg(size(t_bg)-1:size(t_bg)),Pm,tovs(size(tovs))%r,err)
            call polint(x_bg(1,size(t_bg)-1:size(t_bg)),x_bg(3,size(t_bg)-1:size(t_bg)),Pm,tovs(size(tovs))%m,err)
            call polint(x_bg(1,size(t_bg)-1:size(t_bg)),x_bg(4,size(t_bg)-1:size(t_bg)),Pm,tovs(size(tovs))%nu,err)

            deallocate(t_bg,x_bg)
            if (size(pt)-i+1==0) exit
            i=i+1
            r0=tovs(size(tovs))%r
            P0=tovs(size(tovs))%P*(1.d0-eps_Pm)
            m0=tovs(size(tovs))%m
            nu0=tovs(size(tovs))%nu
        enddo

        nuR=0.5d0*log(1.d0-2.d0*G*tovs(size(tovs))%m/tovs(size(tovs))%r/c**2)-tovs(size(tovs))%nu
        tovs(:)%nu=tovs(:)%nu+nuR
        tovs%A=0 !Barotrope

    contains
        function reallocate_bg_a(p,n)
            type(tov_pt_bg), dimension(:), allocatable :: p, reallocate_bg_a
            integer, intent(in) :: n
            integer :: nold,ierr
            allocate(reallocate_bg_a(n),stat=ierr)
            if (ierr /= 0) then
                write(*,*)'reallocate_bg_a: problem in attempt to allocate memory'
                stop
            endif
            if (.not. allocated(p)) return
            nold=size(p,1)
            reallocate_bg_a(1:min(nold,n))=p(1:min(nold,n))
            deallocate(p)
        endfunction reallocate_bg_a
    end subroutine



    subroutine integrate_tov(P0,m0,nu0,r0,Pm,h0,eos_rho)
        use globalVariables,only:G,pi,c
        use odeSolver_tov
        implicit none
        real(8),intent(in)::P0,m0,nu0,r0,Pm
        real(8),intent(inout)::h0
        real(8)::t,x(1:3)
        real(8)::h,emax=1.d-14
        integer:: itmax=100, iflag
        integer::i
        interface
            function eos_rho(p)
                real(8),intent(in)::p
                real(8)::eos_rho
            end function eos_rho
        end interface

        t=r0
        x(1)=P0
        x(2)=m0
        x(3)=nu0
        h=h0

        call rk45ad_ob(tov_eq,tov_cond,t,x,h,itmax,emax,iflag)

        if (allocated(t_bg)) deallocate(t_bg)
        if (allocated(x_bg)) deallocate(x_bg)
        allocate(t_bg(dataSize),x_bg(4,dataSize))
        t_bg=tp
        x_bg(1,:)=xp(1,:)
        do i=1,dataSize
            x_bg(2,i)=eos_rho(xp(1,i))
        enddo
        x_bg(3:4,:)=xp(2:3,:)
        deallocate(tp,xp)

        h0=h

    contains
        subroutine tov_eq(t,x,f)
            real(8), intent(in):: t, x(:)
            real(8), intent(out):: f(:)
            real(8)::r,rho,p,m,nu
            r=t;p=x(1);m=x(2);nu=x(3);rho=eos_rho(p)
            f(1)=-G/r**2*(rho+p/c**2)*(m+4.d0*pi*r**3*P/c**2)/(1.d0-2.d0*G*m/c**2/r)
            f(2)=4.d0*pi*rho*r**2
            f(3)=G/c**2/r**2*(m+4.d0*pi*r**3*P/c**2)/(1.d0-2.d0*G*m/c**2/r)
        endsubroutine tov_eq
        function tov_cond(t,x)
            real(8), intent(in):: t, x(:)
            logical:: tov_cond
            tov_cond = (x(1)<=Pm)
        endfunction tov_cond
    end subroutine integrate_tov

    subroutine tov_pt_get_mr(rhoi,rhof,nTot,eps,eos_pt,eos_p,eos_rho,eos_ga_r,rho,m,r)
!        use io_mod,only:reallocate_a
        real(8),intent(in)::rhoi,rhof,eps
        integer,intent(in)::nTot
        real(8),allocatable,intent(out)::rho(:),m(:),r(:)
        type(tov_pt_bg),allocatable:: tovs(:)
        real(8),allocatable::pt(:)
        real(8)::drho
        integer,allocatable::pt_i(:)
        integer::i
        interface
            function eos_pt(rhoc)
                real(8),intent(in)::rhoc
                real(8),allocatable::eos_pt(:)
            end function
            function eos_rho(p)
                real(8),intent(in)::p
                real(8)::eos_rho
            end function eos_rho
            function eos_p(rho)
                real(8),intent(in)::rho
                real(8)::eos_p
            end function eos_p
            function eos_ga_r(p)
                real(8),intent(in)::p
                real(8)::eos_ga_r
            end function eos_ga_r
        end interface
!        allocate(rho(20),m(20),r(20))
        allocate(rho(nTot),m(nTot),r(nTot))
        drho=(rhof-rhoi)/nTot
        do i=1,nTot
            rho(i)=rhoi+drho*i
            pt=eos_pt(rho(i))
            call tov_pt_driver(rho(i),pt,eps,eos_p,eos_rho,eos_ga_r,tovs,pt_i)
            m(i)=tovs(size(tovs))%m
            r(i)=tovs(size(tovs))%r
        end do
    end subroutine

    function tov_pt_get_rhoc(rho1,rho2,mTarget,pt,eps,eos_p,eos_rho,eos_ga_r)
        use nr_mod,only:rtsec
        real(8),intent(in)::rho1,rho2,mTarget,pt(:),eps
        real(8)::tov_pt_get_rhoc
        real(8)::xacc=1.d-4
        interface
            function eos_rho(p)
                real(8),intent(in)::p
                real(8)::eos_rho
            end function eos_rho
            function eos_p(rho)
                real(8),intent(in)::rho
                real(8)::eos_p
            end function eos_p
            function eos_ga_r(p)
                real(8),intent(in)::p
                real(8)::eos_ga_r
            end function eos_ga_r
        end interface
        tov_pt_get_rhoc=rtsec(m_rho,rho1,rho2,xacc)
    contains
        function m_rho(rho)
            implicit none
            real(8), intent(in) :: rho
            type(tov_pt_bg),allocatable:: tovs(:)
            integer,allocatable:: pt_i(:)
            real(8) :: m_rho
            call tov_pt_driver(rho,pt,eps,eos_p,eos_rho,eos_ga_r,tovs,pt_i)
            m_rho=mTarget-tovs(size(tovs))%m
        end function
    end function
end module

