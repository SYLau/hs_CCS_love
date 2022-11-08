module eos_mod
    !Using EOS table
    !Default format of the EOS table:
    ! - Col 1 and Col 2 are not needed
    ! - Col 3 total energy density in cgs unit
    ! - Col 4 pressure in cgs unit
    implicit none
    private::lnint
    public::eostab
    public::mutab
    public::eos_read,eos_pt,eos_p,eos_rho,eos_ga_r
    public::eos_ga_r_sm
    public::mu_read,eos_mu,eos_chic

    real(8),allocatable::eostab(:,:) !First entry: Column; Second entry: Row
    real(8),allocatable::mutab(:,:) !First entry: Column; Second entry: Row

    !Setting the format of the EOS table
    integer,parameter::crho=3 !Column number of rho (total energy density in cgs unit)
    integer,parameter::cP=4!Column number of P (pressure in cgs unit)
    contains

    subroutine eos_read(path)
    use io_mod,only:file_constTable
        character(len=*)::path
        call file_constTable(path,eostab)
    end subroutine

    function eos_pt(rhoc,errn) !Look for 1st order phase transitions in the EOS table
        real(8),intent(in)::rhoc
        real(8),allocatable::eos_pt(:),temp(:)
        integer,optional,intent(out)::errn
        integer::i,npt
        allocate(temp(size(eostab,2)))
        npt=0
        do i=1,size(eostab,2)-1
            if (eostab(cP,i)==eostab(cP,i+1).and.rhoc>=eostab(crho,i+1)) then
                npt=npt+1
                temp(npt)=eostab(cP,i)
            endif
            if (eostab(cP,i)==eostab(cP,i+1).and.rhoc>eostab(crho,i).and.rhoc<eostab(crho,i+1)) then
                print*,'err: eos_mod eos_pt: rhoc lies within density gap'
                if (present(errn)) then
                    allocate(eos_pt(0))
                    errn=1
                    return
                else
                    stop
                end if
            endif
        end do
        allocate(eos_pt(npt))
        if (npt/=0) eos_pt(1:npt)=temp(npt:1:-1) !reverse order
        if (present(errn)) errn=0
    end function

    function eos_p(rho)
        real(8),intent(in)::rho
        real(8)::eos_p
        integer::i
        do i=1,size(eostab,2)-1
            if (eostab(crho,i)<rho .and. eostab(crho,i+1)>=rho) exit
        enddo
        if (i==size(eostab,2)) then
            print*,'err:eos_p'
            stop
        endif
        eos_p= lnint(eostab(crho,i:i+1),eostab(cP,i:i+1),rho)
    end function eos_p

    function eos_rho(p)
        real(8),intent(in)::p
        real(8)::eos_rho
        integer::i
        do i=1,size(eostab,2)-1
            if (eostab(cP,i)<p .and. eostab(cP,i+1)>=p) exit
        enddo
        if (i==size(eostab,2)) then
            print*,'err:eos_rho'
            stop
        endif
        eos_rho= lnint(eostab(cP,i:i+1),eostab(crho,i:i+1),p)
    end function eos_rho

    function eos_ga_r(p)
    use globalVariables,only:c
        real(8),intent(in)::p
        real(8)::eos_ga_r
        integer::i
        do i=1,size(eostab,2)-1
            if (eostab(cP,i)<p .and. eostab(cP,i+1)>=p) exit
        enddo
        if (i==size(eostab,2)) then
            print*,'err:eos_ga_r'
            stop
        endif
        eos_ga_r= (log(eostab(cP,i+1))-log(eostab(cP,i))) &
        /(log(eostab(crho,i+1))-log(eostab(crho,i)))*(1.d0+p/eos_rho(p)/c**2)
    end function eos_ga_r

    function eos_ga_r_sm(p) !use logistic function to smoothen gamma of the EOS table
    use globalVariables,only:c
        real(8),intent(in)::p
        real(8)::eos_ga_r_sm
        real(8),allocatable,dimension(:)::pm
        real(8)::p1,p2,p3
        real(8)::rho1,rho2,rho3
        real(8)::ga1,ga2,ga_p
        real(8)::width
        integer::i
        allocate(pm(size(eostab,2)-1))
        do i=1,size(eostab,2)-1
            pm(i)=(eostab(cP,i)+eostab(cP,i+1))/2
        end do
        do i=1,size(eostab,2)-2
            if (pm(i)<p .and. pm(i+1)>=p) exit
        enddo
        if (i==size(eostab,2)) then
            print*,'err: eos_ga_r_sm'
            stop
        endif
        p1=eostab(cP,i);p2=eostab(cP,i+1);p3=eostab(cP,i+2)
        rho1=eostab(crho,i);rho2=eostab(crho,i+1);rho3=eostab(crho,i+2)
        ga1=(log(p2)-log(p1))/(log(rho2)-log(rho1))
        ga2=(log(p3)-log(p2))/(log(rho3)-log(rho2))
        if (p2==p3) then
            ga_p=ga1
        elseif (p1==p2) then
            ga_p=ga2
        else
!            width=0.01d0*(p3-p1)
            width=0.005d0*(p3-p1)
            ga_p=ga1+ (ga2-ga1)/(1.d0+exp(-(p-p2)/width)) !logistic function: 1/(1+e^{-(x-x0)/width})
        end if
        eos_ga_r_sm=ga_p*(1.d0+p/eos_rho(p)/c**2)
    end function eos_ga_r_sm

    !Interpolate shear modulus

    subroutine mu_read(path)
    use io_mod,only:file_constTable
        character(len=*)::path
        call file_constTable(path,mutab)
    end subroutine

    impure elemental function eos_mu(p)
        real(8),intent(in)::p
        real(8)::eos_mu
        integer::i
        do i=1,size(mutab,2)-1
            if (mutab(1,i)<p .and. mutab(1,i+1)>=p) exit
        enddo
        if (i==size(mutab,2)) then
            print*,'err:eos_mu'
            stop
        endif
        if (mutab(2,i)==0 .and.mutab(2,i+1)==0) then
            eos_mu=0
        elseif (mutab(2,i)==0) then
            eos_mu=mutab(2,i+1)
        else
            eos_mu=lnint(mutab(1,i:i+1),mutab(2,i:i+1),p)
        end if
    end function eos_mu

    impure elemental function eos_chic(p)
        real(8),intent(in)::p
        real(8)::eos_chic
        integer::i
        do i=1,size(mutab,2)-1
            if (mutab(1,i)<p .and. mutab(1,i+1)>=p) exit
        enddo
        if (i==size(mutab,2)) then
            print*,'err:eos_chic'
            stop
        endif
        if (mutab(3,i)==0 .and.mutab(3,i+1)==0) then
            eos_chic=0
        elseif (mutab(3,i)==0) then
            eos_chic=mutab(3,i+1)
        else
            eos_chic=lnint(mutab(1,i:i+1),mutab(3,i:i+1),p)
        end if
    end function eos_chic

    function lnint(xa,ya,x) !linear log interpolation
        real(8),intent(in)::xa(2),ya(2),x
        real(8)::lnint
        lnint=exp((log(ya(2))-log(ya(1)))/(log(xa(2))-log(xa(1)))*(log(x)-log(xa(1)))+log(ya(1)))
    end function lnint
end module
