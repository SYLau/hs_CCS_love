module odeSolver_tov
    implicit none
    real(8), dimension(:), allocatable :: tp
    real(8), dimension(:,:), allocatable :: xp
    real(8):: tSAS
    integer:: iStep, dataSize

contains

    subroutine rk45ad_ob(fcn,cond,t,x,h,itmax,emax,iflag,xScaleIn,dtSaveIn,hmin) !Solving open boundary ODEs
        implicit none
        integer, intent(in):: itmax
        real(8),intent(in)::emax
        real(8), intent(inout):: t,x(:),h
        integer, intent(out):: iflag
        real(8), optional, intent(in) :: xScaleIn(:)
        real(8), optional, intent(in):: dtSaveIn,hmin
!        real(8) ::delta=0.5e-14
        real(8)::xsave(1:size(x)),tsave, e, xScale(1:size(x)), f(1:size(x)), dtSave
        integer :: n, k, i
        real(8), parameter:: safety=0.9, pGrow= -0.2, pShrink= -0.25, errcon=1.89e-4, tiny= 1.e-30 !The value errcon equals (5/safety)**(1/pGrow), requires scale-up factor to be at most 5

        interface deriv
            subroutine fcn(t,x,f)
                real(8), intent(in):: t, x(:)
                real(8), intent(out):: f(:)
            endsubroutine fcn
        endinterface deriv
        interface condition
            function cond(t,x)
                real(8), intent(in):: t, x(:)
                logical:: cond
            endfunction cond
        endinterface
        n = size(x)
        iflag = 1 ! iflag 0  :-> integration finished
        k = 0
        iStep = 0
        dataSize = 0
        dtSave = 0.d0
        if (present(dtSaveIn)) dtSave = dtSaveIn
        tSAS = t-2.d0*dtSave
        if (allocated(tp)) deallocate(tp)
        if (allocated(xp)) deallocate(xp)
        allocate(tp(256))
        allocate(xp(size(x),size(tp)))

        do
            k = k + 1
            if (k > itmax) exit
            if (cond(t,x)) then !boundary condition matched; exit and save last step
                iflag = 0 !signals integration finished
                call truncate_data
                exit
            end if
            xsave = x
            tsave = t

            if (present(xScaleIn)) then
                xScale(1:n) = xScaleIn(1:n)
            else
                call fcn(t, x, f)
                xScale(1:n) = abs(x(1:n))+abs(h*f(1:n)) + tiny !modify this line to give a custom xScale
            endif

            call rkck45(fcn,t,x,h,xScale(1:n),e) !integrate a single step using Runge-Kutta-Cash-Karp method

            if (e > emax) then !integrate again with shrinked step size
                h= sign(max(abs(safety*h*(e/emax)**pShrink), abs(h)/10.d0), h)
                x = xsave
                t = tsave
                if(t == t+h) pause 'stepsize underflow rk45ad'
            else
                iStep = iStep +1
                k = 0
                if (abs(t-tSAS) > abs(dtSave)) call save_a_step
                if (e/emax > errcon) then !check if the scale-up process is too fast
                    h= safety*h*(e/emax)**pGrow !next step size is scaled up
                else
                    h= 5.d0*h !next step size is scaled up by at most 5
                endif
            endif
            if (present(hmin)) then
                if (abs(h)<abs(hmin)) pause 'err:rk45ad stepsize smaller than hmin'
            endif
        enddo

    contains
        subroutine save_a_step
            use io_mod,only:reallocate_a
            dataSize = dataSize +1
            if (size(tp)<dataSize) then
                tp=reallocate_a(tp,2*size(tp))
                xp=reallocate_a(xp,size(xp,1),2*size(tp))
            endif
            tp(dataSize)=t
            xp(:,dataSize)=x(:)
            tSAS = t
        end subroutine
        subroutine truncate_data
            use io_mod,only:reallocate_a
            tp=reallocate_a(tp,dataSize)
            xp=reallocate_a(xp,size(xp,1),dataSize)
        end subroutine
    end subroutine rk45ad_ob

    subroutine rkck45(fcn,t,x,h,xScale,e)   !Runge-Kutta-Cash-Karp method
        implicit none
        real(8), intent(inout):: t,h,e,x(:),xScale(:)
        integer :: n, i
        real(8), dimension(1:size(x)):: xi, f, f1,f2,f3,f4,f5,f6, x5, err
        real(8):: c20,c21,c30,c31,c32,c40,c41,c42,c43,c51,c52,c53,c54,c60,c61,c62,c63,c64,c65,a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6
        interface deriv
            subroutine fcn(t,x,f)
                real(8), intent(in):: t, x(:)
                real(8), intent(out):: f(:)
            endsubroutine fcn
        endinterface deriv

        c20=0.2d0; c21=0.2d0; c30=0.3d0; c31=0.075d0; c32=0.225d0
        c40=3.d0/5.d0; c41=3.d0/10.d0; c42=-9.d0/10.d0; c43=6.d0/5.d0
        c51=-11.d0/54.d0; c52=5.d0/2.d0; c53=-70.d0/27.d0; c54=35.d0/27.d0
        c60=0.875d0; c61=1631.d0/55296.d0; c62=175.d0/512.d0; c63=575.d0/13824.d0; c64=44275.d0/110592.d0; c65=253.d0/4096.d0
        a1=2825.d0/27648.d0; a2=0.d0; a3=18575.d0/48384.d0; a4=13525.d0/55296.d0; a5=277.d0/14336.d0; a6=0.25d0
        b1=37.d0/378.d0; b2=0.d0; b3=250.d0/621.d0; b4=125.d0/594.d0; b5=0.d0; b6=512.d0/1771.d0
        n = size(x)
        xi(1:n) = x(1:n)
        call fcn(t, xi, f)
        f1(1:n) = h*f(1:n)

        xi(1:n) = x(1:n) + c21*f1(1:n)
        call fcn(t+ c20*h, xi, f)
        f2(1:n) = h*f(1:n)

        xi(1:n) = x(1:n) + c31*f1(1:n) + c32*f2(1:n)
        call fcn(t+ c30*h, xi, f)
        f3(1:n) = h*f(1:n)

        xi(1:n) = x(1:n) + c41*f1(1:n) + c42*f2(1:n) + c43*f3(1:n)
        call fcn(t+ c40*h, xi, f)
        f4(1:n) = h*f(1:n)

        xi(1:n) = x(1:n) + c51*f1(1:n) + c52*f2(1:n) + c53*f3(1:n) + c54*f4(1:n)
        call fcn(t+h, xi, f)
        f5(1:n) = h*f(1:n)

        xi(1:n) = x(1:n) + c61*f1(1:n) + c62*f2(1:n) + c63*f3(1:n) + c64*f4(1:n) + c65*f5(1:n)
        call fcn(t+ c60*h, xi, f)
        f6(1:n) = h*f(1:n)

        x5(1:n) = x(1:n) + b1*f1(1:n) + b3*f3(1:n) + b4*f4(1:n) + b5*f5(1:n) + b6*f6(1:n)
        x(1:n) = x(1:n) + a1*f1(1:n) + a3*f3(1:n) + a4*f4(1:n) + a5*f5(1:n)  +  a6*f6(1:n)
        t = t + h
        err(1:n) = dabs(x(1:n) - x5(1:n)) /dabs(xScale(1:n))
        e= maxval(err(1:n))
    end subroutine rkck45

endmodule odeSolver_tov

