
program hs_love
    use eos_hs_gen_mod,only:mitb_par,hs_path
    use eos_hs_gen_mod,only:gen_hybrid_star_eos
    implicit none

    type(mitb_par)::mitpa_in
    type(hs_path)::hs_gen_path
    real(8)::pt_qm_nm
    character(len=99)::outdir='./'
    character(len=255)::eos_path,mu_path
    character(len=99)::Gap
    integer::qm_ntot
    real(8)::rho0
    real(8)::mind,maxd, diff, step
    real(8), dimension(100)::lprho0
    integer:: nump, i, j
    integer,allocatable,dimension(:)::ode_steps

    !generate the HS EOS using the modified MIT Bag model (Eq. (3.1) of ): https://arxiv.org/abs/nucl-th/0411016 ; stored in hs_gen_path%hs
    !generate the shear modulus using Eq. (105): https://arxiv.org/abs/hep-ph/0702021 ; stored in hs_gen_path%mu
    mitpa_in%a4=0.85d0; mitpa_in%a2=100.d0**2; mitpa_in%beff=160.d0**4; mitpa_in%gap=25.d0   ! gap 5 - 25 MeV ! shear modulus goes as gap ^2

    !above show the parameters of the shear modulus equation, coefficients a4 -- a2 -- effective bag constant -- gap parameter (energy of LOFF pairing (maybe)

    pt_qm_nm=2.06818d34 !Quark-Hadron phase transition pressure (should be determined using Maxwell's construction)
    !could try higher value of high transition pressure

    qm_ntot=150 !number of lines for the quark matter part of the EOS table (no need to change)


    hs_gen_path%nm='eos_apr.txt'
    !hs_gen_path%nm='eos_SLy4.txt'
    hs_gen_path%hs='eos_hs.txt'
    hs_gen_path%mu='mu_hs.txt'
    !eos-paths
    !
    ! try to code up new model of CSS
    mind=0.5d15
    maxd=4d15
    diff=maxd-mind
    nump=35
    write(*,*)'diff', diff
    step=(maxd-mind)/nump
    write(*,*) 'step', step


    do j=1,nump
        step=(maxd-mind)/nump
        lprho0(j)=mind+(j-1)*step
        !write(*,*)'pho', lprho0(j)
    enddo

    do i =1, nump
    !generates the EOS for hybrid star by combining the nuclear matter (fluid) and quark matter (solid) EOS's
    call gen_hybrid_star_eos(hs_gen_path,mitpa_in,pt_qm_nm,qm_ntot)
    !maybe later on, generate other parameterized equations of state to stitch to fluid and then solve

    !Calculate Love number of a single HS
    eos_path=hs_gen_path%hs !HS EOS (must have at least one phase transition, denoted with two lines with equal p) [see description in cal_single_love_hybrid_star]
    mu_path=hs_gen_path%mu  !CCS core elasticity [see description in cal_single_love_hybrid_star]

    rho0=lprho0(i)
    write(*,*) 'pho', rho0
    !loop through similar central densities


    call cal_single_love_hybrid_star(outdir,rho0,eos_path,mu_path,2.d0,i) !*****Check this subroutine in detail*****

    write(*,*) 'Start Convergence Test'

    !A convergence test to show numerical accuracy:
    !   varies the numerical integration stepsize to see how the result changes
    !   important to keep track of how accurate the solution is
    !Comment it out if not needed
    ode_steps=[750,1500,3000,4500,6000] ! define (half of) number of steps for integration for Love number
    !call run_love_hybrid_star_conv_test(outdir,rho0,eos_path,mu_path,ode_steps,2.d0)
    !using EOS, gamma, using tabulated equation of state, there will be a discrete gamma, uses interpolation
    !derivative for gamma wont be smooth anymore
    !adaptive grid size and uniform grid size
    !adaptive, this can be problematic, but gamma isn't smooth, which can be a problem
    !step size is adaptive to accuracy, the adaptive can't be used anymore because gamma is discrete, use uniform grid integrator
    !now we need to test accuracy by using convergence test

    end do

    contains

    subroutine cal_single_love_hybrid_star(outdir,rhoc,eos_path,mu_path,ell,i)
        !Requirement for the eos_path file (see eos_mod): In cgs units: Col 1: index; Col 2: nb; Col 3: rho; Col 4: p
        !eos_path must have at least one phase transition, denoted with two lines with equal p.
        !Requirement for the mu_path file (see eos_mod): In cgs units: Col 1: p; Col 2: mu; Col 2: chic, where chic = dmu/dp
        !In mu_path, p must be increasing (For first order phase transition, remove line with p=pt at fluid side)
        !only supports ell=2 at this moment
        use globalVariables,only:M_sun
        use globalVariables,only:G,pi
        use globalVariables,only:c
        use eos_mod,only:eos_pt,eos_read,eos_p,eos_rho,eos_ga_r
        use eos_mod,only:eos_ga_r_sm
        use eos_mod,only:mu_read,eos_mu,eos_chic
        use bg_mod,only:load_bgs
        use bg_mod,only:load_bgs_mu
        use io_mod,only:writeArray,append_a
        use gnuPlot,only:gpf_plot2D
        use tov_pt_mod,only:tov_pt_driver,tov_pt_bg
        use love_gr_mod,only:fl_adapt=>adaptive
        use love_gr_mod,only:get_love_gr_fluid
        use love_gr_mod,only:t_lo,x_lo
        use love_gr_mod,only:love_eq_ell2
        use love_gr_s1f2_mod,only:hs_adapt=>adaptive
        use love_gr_s1f2_mod,only:get_love_gr_s1f2
        use love_gr_s1f2_mod,only:t_s,x_s,t_f,x_f
        use ctime_mod,only:print_time
        character(len=*),intent(in)::outdir
        character(len=99)::outfle
        real(8),intent(in)::rhoc !Central density (cgs unit)
        character(len=*),intent(in)::eos_path,mu_path
        real(8),intent(in)::ell !Spherical harmonic degree l

        real(8),allocatable::pt(:) !Stores the pressure of the phase transition points
        integer::err_pt
        type(tov_pt_bg),allocatable:: tovs(:) !Stores the background stellar structure (Solution to TOV equations)
        real(8),allocatable,dimension(:):: mu,chic
        real(8)::mtot,rtot
        integer,allocatable:: pt_i(:) !Stores the position of density discontinuity in the array tovs - phase transition
        real(8)::yout
        real(8)::kl,kl_fluid
        real(8),allocatable,dimension(:)::xplot1,yplot1,xplot2,yplot2
        integer :: i

        call eos_read(eos_path)
        pt=eos_pt(rhoc,err_pt) !******changed
        if (err_pt /= 0) then !err_pt/= 0 means the central density lies within the density gap
            print*,'warning: HS rho0 within density gap at rho0 = ',rhoc
            return
        endif
        if (size(pt)<1) print*,'warning: HS no 1st order phase transition found'

        call print_time(msg='Solving TOV equations')
!        call tov_pt_driver(rhoc,pt,1.d-14,eos_p,eos_rho,eos_ga_r,tovs,pt_i) !Solve TOV equations, store solution in derived type: tovs
        call tov_pt_driver(rhoc,pt,1.d-14,eos_p,eos_rho,eos_ga_r_sm,tovs,pt_i) !smoothen gamma using logistic function

        !Calculate shear modulus
        call mu_read(mu_path)

        if (size(pt)>=1) then !*********change [Nov 04 2022]
            mu=eos_mu(tovs%p)
            chic=eos_chic(tovs%p)
            mu(pt_i(1)+1:size(mu))=0
            chic(pt_i(1)+1:size(chic))=0
        else
            mu=tovs%p*0
            chic=tovs%p*0
        end if


        mtot=tovs(size(tovs))%m
        rtot=tovs(size(tovs))%r

        call load_bgs(tovs%r,tovs%p,tovs%rho,tovs%m,tovs%nu,tovs%ga,tovs%A,pt_i) !load background solution into bg_mod: required for solving the perturbation eqts
        call load_bgs_mu(mu,chic)

        call print_time(msg='Solving for Love number')
        fl_adapt=.false. !switch off adaptive gridsize for integration (unstable when using discrete EOS tables)
        call get_love_gr_fluid(ell,kl_fluid) !Calculate the tidal Love number.
        print*,'Love number (fluid) = ', kl_fluid
        print*,'Tidal deformability (fluid) = ', kl_fluid*2/3*(rtot*c**2/G/mtot)**5


        if (size(pt)>=1) then !******changed
            hs_adapt=.false. !switch off adaptive gridsize for integration (unstable when using discrete EOS tables)
!            print*, 'check 2'
            call get_love_gr_s1f2(ell,.true.,yout) !Calculate the solution of (y_out = R H0'/H0) just outside the surface !*****Check this module in detail*****
!            print*, 'check 1'
            kl=love_eq_ell2(mtot,rtot,yout) !Using y_out (= R H0'/H0) to calculate the tidal Love number of HS with CCS quark core
            print*,'Love number (HS) = ', kl
            print*,'Tidal deformability = ', kl*2/3*(rtot*c**2/G/mtot)**5
        else
            kl=kl_fluid
            print*,'No solid layer in hybrid stars'
        endif

!        print*, 'pass!!!'

        print *, "Output files"
        call gpf_plot2D('density profile','r (cm)','rho (g/cm^3)','set style data lines',tovs%r,tovs%rho & !Generate the density profile plot
            ,filename=trim(outdir)//'rho_r.gp',display=.false.)
        call gpf_plot2D('gamma profile','r (cm)','gamma','set style data lines',tovs%r,tovs%ga & !Generate the adiabatic index profile plot
            ,filename=trim(outdir)//'ga_r.gp',display=.false.)
        call writeArray(trim(outdir)//'bg.dat',[character(len=25)::'r (cm)','P (cgs)','rho (cgs)','m (g)','nu','mu (cgs)'])
        call writeArray(trim(outdir)//'bg.dat',tovs%r,tovs%p,tovs%rho,tovs%m,tovs%nu,mu,st='old',ac='append') !Generate the data file containing the static profile

        outfle='model_info_'//hs_gen_path%nm(1:7)//'.dat'

        if (i==1) then
            call writeArray(trim(outdir)//outfle &
                ,[character(len=25)::'rho0(cgs)','P0(cgs)','M(M_sun)','R(km)','P_min(P0)','rho_min(rho0)' &
                ,'P_qm_nm(cgs)','k_2_fluid','k_2','C'])
        endif

        call writeArray(trim(outdir)//outfle,[tovs(1)%rho],[tovs(1)%P],[mtot/M_sun],[rtot/1.d5] &
            ,[tovs(size(tovs))%P/tovs(1)%P],[tovs(size(tovs))%rho/tovs(1)%rho] &
            ,[pt_qm_nm],[kl_fluid],[kl],[(G/(2.99792458d10**2)*(mtot/rtot))], st='unknown',ac='append') !Generate the data file containing the model information, including Love numbers

        if (size(pt)>=1) then !******changed
            xplot1=t_lo; yplot1=x_lo
            xplot2=append_a(t_s,t_f); yplot2=append_a(t_s,t_f)*append_a(x_s(2,:),x_f(2,:))/append_a(x_s(1,:),x_f(1,:))
            call gpf_plot2D('eta profile','r (cm)','eta','set style data lines',xplot1,yplot1,ls1='title "fluid"' &
                ,x2=xplot2,y2=yplot2,ls2='title "CCS core"',filename=trim(outdir)//'love_eq_eta.gp',display=.false.) !Generate the (r H_0^\prime/H_0) plot
            call writeArray(trim(outdir)//'love_eq_solid_core.dat' &
                ,[character(len=25)::'r (cgs)','H0 (cgs)','H0` (cgs)','K (cgs)','W (cgw)','V (cgw)','T2 (cgs)'])
            call writeArray(trim(outdir)//'love_eq_solid_core.dat',t_s,x_s(1,:),x_s(2,:),x_s(3,:),x_s(4,:),x_s(5,:) &
                ,x_s(6,:),st='old',ac='append') !Generate the data file containing the solid core solution
        endif
    end subroutine

     subroutine run_love_hybrid_star_conv_test(outdir,rhoc,eos_path,mu_path,ode_steps,ell)
        use globalVariables,only:M_sun
        use globalVariables,only:G,pi
        use globalVariables,only:c
        use eos_mod,only:eos_pt,eos_read,eos_p,eos_rho,eos_ga_r
        use eos_mod,only:eos_ga_r_sm
        use eos_mod,only:mu_read,eos_mu,eos_chic
        use bg_mod,only:load_bgs
        use bg_mod,only:load_bgs_mu
        use io_mod,only:writeArray,append_a
        use gnuPlot,only:gpf_plot2D
        use tov_pt_mod,only:tov_pt_driver,tov_pt_bg
        use love_gr_mod,only:fl_adapt=>adaptive,fl_ug_steps=>ug_steps
        use love_gr_mod,only:get_love_gr_fluid
        use love_gr_mod,only:love_eq_ell2
        use love_gr_s1f2_mod,only:hs_adapt=>adaptive,hs_ug_steps=>ug_steps
        use love_gr_s1f2_mod,only:get_love_gr_s1f2
        use ctime_mod,only:print_time
        character(len=*),intent(in)::outdir
        real(8),intent(in)::rhoc !Central density (cgs unit)
        character(len=*),intent(in)::eos_path,mu_path
        integer,dimension(:),intent(in)::ode_steps
        real(8),intent(in)::ell !Spherical harmonic degree l

        real(8),allocatable::pt(:) !Stores the pressure of the phase transition points
        integer::err_pt
        type(tov_pt_bg),allocatable:: tovs(:) !Stores the background stellar structure (Solution to TOV equations)
        real(8),allocatable,dimension(:):: mu,chic
        real(8)::mtot,rtot
        integer,allocatable:: pt_i(:) !Stores the position of density discontinuity in the array tovs
        real(8)::yout
        real(8)::kl,kl_fluid
        real(8),allocatable,dimension(:)::kl_fluid_cnv,kl_cnv
        integer::i


        call eos_read(eos_path)
!        pt=eos_pt(rhoc) !******Original
!        if (size(pt)<1) print*,'warning: HS no 1st order phase transition found' !******Original
        pt=eos_pt(rhoc,err_pt) !******changed
        if (err_pt /= 0) then !err_pt/= 0 means the central density lies within the density gap
            print*,'warning: HS rho0 within density gap at rho0 = ',rhoc
            return
        endif
        if (size(pt)<1) print*,'warning: HS no 1st order phase transition found'

        call print_time(msg='Solving TOV equations')
!        call tov_pt_driver(rhoc,pt,1.d-14,eos_p,eos_rho,eos_ga_r,tovs,pt_i) !Solve TOV equations, store solution in derived type: tovs
        call tov_pt_driver(rhoc,pt,1.d-14,eos_p,eos_rho,eos_ga_r_sm,tovs,pt_i) !smoothen gamma using logistic function

        !Calculate shear modulus
        call mu_read(mu_path)
        if (size(pt)>=1) then !*********change [Nov 04 2022]
            mu=eos_mu(tovs%p)
            chic=eos_chic(tovs%p)
            mu(pt_i(1)+1:size(mu))=0
            chic(pt_i(1)+1:size(chic))=0
        else
            mu=tovs%p*0
            chic=tovs%p*0
        end if

        mtot=tovs(size(tovs))%m
        rtot=tovs(size(tovs))%r

        call load_bgs(tovs%r,tovs%p,tovs%rho,tovs%m,tovs%nu,tovs%ga,tovs%A,pt_i) !load background solution into bg_mod
        call load_bgs_mu(mu,chic)

        allocate(kl_fluid_cnv(0),kl_cnv(0))
        fl_adapt=.false.
        hs_adapt=.false.
        do i=1,size(ode_steps)
            call print_time(msg='Convergence test: Solving for Love number')
            fl_ug_steps=ode_steps(i) ! set (half of) number of steps in integration
            hs_ug_steps=ode_steps(i)
            call get_love_gr_fluid(ell,kl_fluid) !Calculate the tidal Love number.
!            if (eos_p(rhoc) >= tovs(pt_i(1))%p) then !******Original
            if (size(pt)>=1) then !******changed
                call get_love_gr_s1f2(ell,.true.,yout) !Calculate the tidal Love number of HS with CCS quark core.
                kl=love_eq_ell2(mtot,rtot,yout)
            else
                kl=kl_fluid
            endif

            kl_fluid_cnv=append_a(kl_fluid_cnv,[kl_fluid])
            kl_cnv=append_a(kl_cnv,[kl])
        enddo

        kl_fluid_cnv=(kl_fluid_cnv-kl_fluid_cnv(size(kl_fluid_cnv)))/kl_fluid_cnv !calculate fractional change
        kl_cnv=(kl_cnv-kl_cnv(size(kl_cnv)))/kl_cnv

        print *, "Output files"
        call gpf_plot2D('Love number convergence','steps','abs frac change' &
            ,options='set logscale xy' &
            ,x1=real(ode_steps,kind=8),y1=abs(kl_fluid_cnv),ls1='title "Fluid"' &
            ,x2=real(ode_steps,kind=8),y2=abs(kl_cnv),ls2='title "CCS"' &
            ,filename=trim(outdir)//'kl_conv.gp',display=.false.)

        call writeArray(trim(outdir)//'model_info.dat' &
            ,[character(len=25)::'rho0 (cgs)','P0 (cgs)','M (M_sun)','R (km)','P_min (P0)','rho_min (rho0)' &
            ,'P_qm_nm (cgs)','k_2_fluid','k_2'])
        call writeArray(trim(outdir)//'model_info.dat',[tovs(1)%rho],[tovs(1)%P],[mtot/M_sun],[rtot/1.d5] &
            ,[tovs(size(tovs))%P/tovs(1)%P],[tovs(size(tovs))%rho/tovs(1)%rho] &
            ,[pt_qm_nm],[kl_fluid],[kl],st='old',ac='append') !Generate the data file containing the model information, including Love numbers

    end subroutine

end program

