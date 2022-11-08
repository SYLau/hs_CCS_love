module eos_hs_gen_mod
    private
    public::mitb_par
    public::hs_path
    public::gen_hybrid_star_eos

    type mitb_par
        real(8)::a4,a2,beff,gap
    end type mitb_par

    type hs_path
        character(len=255)::nm,hs,mu
    end type hs_path

    contains

    subroutine gen_hybrid_star_eos(inpath,mitpa_in,pt,qm_ntot) !Generate the hybrid star EOS by connecting NM EOS and QM EOS (modified Bag model); assuming 1st order phase transition
        use eos_mod,only:eos_read,eostab
        use eos_mitb_mod,only:mitpa,eos_rho,eos_nq
        use eos_mitb_mod,only:eos_mu,eos_chic
        type(hs_path),intent(in)::inpath
        type(mitb_par),intent(in)::mitpa_in
        real(8),intent(in)::pt
        integer,intent(in)::qm_ntot
        integer,parameter::cindx=1,cnb=2,crho=3,cP=4
        integer,allocatable,dimension(:)::indx
        real(8),allocatable,dimension(:)::rho,p,nb
        real(8)::pmax
        integer::i
        integer::nm_ntot,ntot
        integer::uni

        call eos_read(inpath%nm)
        do i=1,size(eostab,2)
            if (eostab(cP,i)>=pt) exit
        end do
        if (i==size(eostab,2)) then
            write(*,*) 'err: gen_hybrid_star_eos Pt higher than NM EOS'
            stop
        end if

        !pt is phase transition
        !look in eos table when the pressure is greater than or equal to the transition pressure

        nm_ntot=i
        ntot=qm_ntot+nm_ntot
        allocate(indx(ntot),nb(ntot),rho(ntot),p(ntot))
        indx(1:nm_ntot-1)=nint(eostab(cindx,1:nm_ntot-1))
        nb(1:nm_ntot-1)=eostab(cnb,1:nm_ntot-1)
        rho(1:nm_ntot-1)=eostab(crho,1:nm_ntot-1)
        p(1:nm_ntot-1)=eostab(cp,1:nm_ntot-1)

        indx(nm_ntot)=indx(nm_ntot-1)+1
        p(nm_ntot)=pt
        rho(nm_ntot)=lnint(eostab(cp,nm_ntot-1:nm_ntot),eostab(crho,nm_ntot-1:nm_ntot),pt) !loglin interpolation
        nb(nm_ntot)=lnint(eostab(cp,nm_ntot-1:nm_ntot),eostab(cnb,nm_ntot-1:nm_ntot),pt)

        mitpa%a4=mitpa_in%a4; mitpa%a2=mitpa_in%a2; mitpa%beff=mitpa_in%beff; mitpa%gap=mitpa_in%gap
        !renaming paramters a4, a2, mit bag constant and gap parameter
        indx(nm_ntot+1:ntot)=2*indx(nm_ntot)
        p(nm_ntot+1)=pt
        rho(nm_ntot+1)=eos_rho(pt)
        nb(nm_ntot+1)=eos_nq(pt)/3
        pmax=eostab(cp,size(eostab,2))

        do i=2,qm_ntot
            p(nm_ntot+i)=exp((log(pmax)-log(pt))/qm_ntot*i+log(pt))
            rho(nm_ntot+i)=eos_rho(p(nm_ntot+i))
            nb(nm_ntot+i)=eos_nq(p(nm_ntot+i))/3
        end do

        open(newunit=uni,file=inpath%hs,status='replace')
        do i=1,ntot
            write(uni,'(i10,3(es25.16))')indx(i),nb(i),rho(i),p(i)
        end do
        close(uni)

        open(newunit=uni,file=inpath%mu,status='replace')
        do i=1,nm_ntot-1 !skipping i = nm_ntot, since p(nm_ntot)=p(nm_ntot+1)
            write(uni,'(3(es25.16))')p(i),0.d0,0.d0
        end do
        do i=nm_ntot+1,ntot
            write(uni,'(3(es25.16))')p(i),eos_mu(p(i)),eos_chic(p(i))
        end do
        close(uni)
    end subroutine
    function lnint(xa,ya,x) !linear log interpolation
        real(8),intent(in)::xa(2),ya(2),x
        real(8)::lnint
        lnint=exp((log(ya(2))-log(ya(1)))/(log(xa(2))-log(xa(1)))*(log(x)-log(xa(1)))+log(ya(1)))
    end function lnint
end module eos_hs_gen_mod
