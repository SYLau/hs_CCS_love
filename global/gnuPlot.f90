module gnuPlot
    !Version 2.0 - Interface with the object oriented GnuPlot code (ogpf)
    implicit none
    private
    public::gpf_plot2D
    public::gpf_surf

    contains

    subroutine gpf_plot2D(title,xlabel,ylabel,options,x1,y1,ls1,x2,y2,ls2,x3,y3,ls3 &
    ,x4,y4,ls4,filename,display)
        use ogpf, only:wp,gpf
        type(gpf)::gp
        character(len=*),intent(in)::title,xlabel,ylabel,options
        real(wp),dimension(:),intent(in)::x1,y1
        character(len=*),intent(in),optional::ls1
        real(wp),dimension(:),intent(in),optional::x2,y2
        character(len=*),intent(in),optional::ls2
        real(wp),dimension(:),intent(in),optional::x3,y3
        character(len=*),intent(in),optional::ls3
        real(wp),dimension(:),intent(in),optional::x4,y4
        character(len=*),intent(in),optional::ls4
        character(len=*),intent(in),optional::filename
        logical,intent(in),optional::display
!        logical::label1=.false.,data2=.false.,label2=.false.,data3=.false.,label3=.false. &
!        ,data4=.false.,label4=.false.
        call gp%reset()
        if (present(filename))call gp%filename(filename)
        if (present(display))call gp%display(display)
        call gp%title(title)
        call gp%xlabel(xlabel)
        call gp%ylabel(ylabel)
        call gp%options(options)

        call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2,ls2=ls2,x3=x3,y3=y3,ls3=ls3 &
        ,x4=x4,y4=y4,ls4=ls4)

!        if (present(ls1)) label1=.true.
!        if (present(x2).and.present(y2)) data2=.true.
!        if (present(ls2)) label2=.true.
!        if (present(x3).and.present(y3)) data3=.true.
!        if (present(ls3)) label3=.true.
!        if (present(x4).and.present(y4)) data4=.true.
!        if (present(ls4)) label4=.true.
!        if (.not.(data2.or.data3.or.data4)) then
!            if (.not.label1) then
!                call gp%plot(x1,y1)
!            else
!                call gp%plot(x1,y1,ls1=ls1)
!            end if
!        elseif (.not.(data3.or.data4)) then
!            if (.not.(label1.or.label2)) then
!                call gp%plot(x1,y1,x2=x2,y2=y2)
!            elseif (.not.label2) then
!                call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2)
!            elseif (.not.label1) then
!                call gp%plot(x1,y1,x2=x2,y2=y2,ls2=ls2)
!            else
!                call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2,ls2=ls2)
!            end if
!        elseif (.not.(data4)) then
!            if (.not.(label1.or.label2.or.label3)) then
!                call gp%plot(x1,y1,x2=x2,y2=y2,x3=x3,y3=y3)
!            elseif (.not.(label2.or.label3)) then
!                call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2,x3=x3,y3=y3)
!            elseif (.not.(label3.or.label1)) then
!                call gp%plot(x1,y1,x2=x2,y2=y2,ls2=ls2,x3=x3,y3=y3)
!            elseif (.not.(label1.or.label2)) then
!                call gp%plot(x1,y1,x2=x2,y2=y2,x3=x3,y3=y3,ls3=ls3)
!            elseif (.not.label3) then
!                call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2,ls2=ls2,x3=x3,y3=y3)
!            elseif (.not.label1) then
!                call gp%plot(x1,y1,x2=x2,y2=y2,ls2=ls2,x3=x3,y3=y3,ls3=ls3)
!            elseif (.not.label2) then
!                call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2,x3=x3,y3=y3,ls3=ls3)
!            else
!                call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2,ls2=ls2,x3=x3,y3=y3,ls3=ls3)
!            end if
!        else
!            if (.not.(label1.or.label2.or.label3.or.label4)) then
!                call gp%plot(x1,y1,x2=x2,y2=y2,x3=x3,y3=y3,x4=x4,y4=y4)
!            elseif (.not.(label2.or.label3.or.label4)) then
!                call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2,x3=x3,y3=y3,x4=x4,y4=y4)
!            elseif (.not.(label1.or.label3.or.label4)) then
!                call gp%plot(x1,y1,x2=x2,y2=y2,ls2=ls2,x3=x3,y3=y3,x4=x4,y4=y4)
!            elseif (.not.(label1.or.label2.or.label4)) then
!                call gp%plot(x1,y1,x2=x2,y2=y2,x3=x3,y3=y3,ls3=ls3,x4=x4,y4=y4)
!            elseif (.not.(label1.or.label2.or.label3)) then
!                call gp%plot(x1,y1,x2=x2,y2=y2,x3=x3,y3=y3,x4=x4,y4=y4,ls4=ls4)
!            elseif (.not.(label3.or.label4)) then
!                call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2,ls2=ls2,x3=x3,y3=y3,x4=x4,y4=y4)
!            elseif (.not.(label2.or.label4)) then
!                call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2,x3=x3,y3=y3,ls3=ls3,x4=x4,y4=y4)
!            elseif (.not.(label2.or.label3)) then
!                call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2,x3=x3,y3=y3,x4=x4,y4=y4,ls4=ls4)
!            elseif (.not.(label1.or.label4)) then
!                call gp%plot(x1,y1,x2=x2,y2=y2,ls2=ls2,x3=x3,y3=y3,ls3=ls3,x4=x4,y4=y4)
!            elseif (.not.(label1.or.label3)) then
!                call gp%plot(x1,y1,x2=x2,y2=y2,ls2=ls2,x3=x3,y3=y3,x4=x4,y4=y4,ls4=ls4)
!            elseif (.not.(label1.or.label2)) then
!                call gp%plot(x1,y1,x2=x2,y2=y2,x3=x3,y3=y3,ls3=ls3,x4=x4,y4=y4,ls4=ls4)
!            elseif (.not.label1) then
!                call gp%plot(x1,y1,x2=x2,y2=y2,ls2=ls2,x3=x3,y3=y3,ls3=ls3,x4=x4,y4=y4 &
!                ,ls4=ls4)
!            elseif (.not.label2) then
!                call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2,x3=x3,y3=y3,ls3=ls3,x4=x4,y4=y4 &
!                ,ls4=ls4)
!            elseif (.not.label3) then
!                call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2,ls2=ls2,x3=x3,y3=y3,x4=x4,y4=y4 &
!                ,ls4=ls4)
!            elseif (.not.label4) then
!                call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2,ls2=ls2,x3=x3,y3=y3,ls3=ls3 &
!                ,x4=x4,y4=y4)
!            else
!                call gp%plot(x1,y1,ls1=ls1,x2=x2,y2=y2,ls2=ls2,x3=x3,y3=y3,ls3=ls3 &
!                ,x4=x4,y4=y4,ls4=ls4)
!            end if
!        end if

    end subroutine

    subroutine gpf_surf(title,xlabel,ylabel,options,x1,y1,z1,palette,filename,display)
        use ogpf, only:wp,gpf
        type(gpf)::gp
        character(len=*),intent(in)::title,xlabel,ylabel,options
        real(wp),dimension(:,:),intent(in)::x1,y1,z1
        character(len=*),intent(in),optional::palette
        character(len=*),intent(in),optional::filename
        logical,intent(in),optional::display

        call gp%reset()
        if (present(filename))call gp%filename(filename)
        if (present(display))call gp%display(display)
        call gp%title(title)
        call gp%xlabel(xlabel)
        call gp%ylabel(ylabel)
        call gp%options(options)
!        call gp%options('unset key')
!        call gp%options('set pm3d map')
!        call gp%options('set logscale zcb')
!        call gp%options('set cblabel "---"')
!        call gp%options('set rmargin at screen 0.8')
!        call gp%options('set lmargin at screen 0.1')

        call gp%surf(x1,y1,z1,palette)

    end subroutine gpf_surf

end module gnuPlot

!module gnuPlot_v1_0
!    !Version 1.0 - Self written interface to plot with gnuPlot. Has limited functionality.
!    !
!    ! Sample gnu command to plot with premade plt files:
!    !          call system('cd plots && gnuplot -persist "plotOrbitKeplerAna.plt"  ')
!    ! gnu command:
!    !        set title "Analytic Kepler Orbit"
!    !        set xlabel "x / AU"
!    !        set ylabel "y / AU"
!    !        plot '../output/keplerOrbitAna.txt' using 1:2 with lines title "Keplerian orbit"
!    implicit none
!
!    type str
!        character(len=255)::d
!    end type
!
!    interface gnuPlot2D
!        module procedure gnuPlot2D_fx, gnuPlot2D_xy, gnuPlot2D_xy2, gnuPlot2D_xy_list
!    endinterface gnuPlot2D
!
!contains
!
!    !Plot y = f(x) graphs
!    subroutine gnuPlot2D_fx(fcn, Dir, title, xlabel, ylabel, xi, xf, dataSize, display,output,plot_title_in)
!        real(8), intent(in) :: xi, xf
!        integer, intent(in):: dataSize
!        character (len=*) :: Dir, title, xlabel, ylabel
!        character (len=*),optional :: plot_title_in
!        character (len=255):: plot_title,command
!        logical:: display,output
!        type(str)::title_str(1)
!        real(8) :: x(dataSize), y(dataSize), d
!        integer :: i
!        interface
!            function fcn(x)
!                implicit none
!                real(8), intent(in) :: x
!                real(8) :: fcn
!            endfunction fcn
!        endinterface
!        title_str(1)%d=title
!        plot_title=title
!        if (present(plot_title_in)) plot_title=plot_title_in
!        d = (xf-xi)/dataSize
!        do i = 1, dataSize
!            x(i) = xi + d*i
!            y(i) = fcn(x(i))
!        enddo
!        call gnuPlot2D_xy(Dir, title, xlabel, ylabel, x,y,display,output,plot_title_in=plot_title_in)
!    endsubroutine gnuPlot2D_fx
!
!    !Plot x-y graphs
!    subroutine gnuPlot2D_xy(Dir, title, xlabel, ylabel, x,y,display,output,plot_title_in)
!        real(8), intent(in), dimension(:) :: x,y
!        character (len=*) :: Dir, title, xlabel, ylabel
!        character (len=255):: plot_title,command
!        character (len=*),optional :: plot_title_in
!        logical:: display,output
!        type(str)::title_str(1)
!        real(8)::y_tab(size(y),1)
!        title_str(1)%d=title
!        y_tab(:,1)=y
!        plot_title=title
!        if (present(plot_title_in)) plot_title=plot_title_in
!        call gnuCommandFile(Dir, plot_title,title_str, xlabel, ylabel, display,output)
!        call gnuDataFile(Dir, x,y_tab)
!        write(command,*)'cd '//trim(Dir)//' && gnuplot -persist gnuCF_f90.plt'//' &'
!        call execute_command_line(trim ( command ) )
!        call gnuPlot_clean(Dir)
!    endsubroutine gnuPlot2D_xy
!
!    !Plot x-yy graphs
!    subroutine gnuPlot2D_xy2(Dir, title1, title2, xlabel, ylabel, x,y1,y2,display,output,plot_title_in)
!        real(8), intent(in), dimension(:) :: x,y1,y2
!        character (len=*) :: Dir, title1, title2, xlabel, ylabel
!        character (len=255):: plot_title,command
!        character (len=*),optional :: plot_title_in
!        logical:: display,output
!        type(str)::title_str(2)
!        real(8)::y_tab(size(y1),2)
!        title_str(1)%d=title1;title_str(2)%d=title2
!        y_tab(:,1)=y1;y_tab(:,2)=y2
!        plot_title=title1
!        if (present(plot_title_in)) plot_title=plot_title_in
!        call gnuCommandFile(Dir, plot_title,title_str, xlabel, ylabel, display,output)
!        call gnuDataFile(Dir, x,y_tab)
!        write(command,*)'cd '//trim(Dir)//' && gnuplot -persist gnuCF_f90.plt'//' &'
!        call execute_command_line(trim ( command ) )
!        call gnuPlot_clean(Dir)
!    endsubroutine gnuPlot2D_xy2
!
!    subroutine gnuPlot2D_xy_list(Dir, title, xlabel, ylabel, x,y,display,output,plot_title_in)
!        real(8), intent(in):: x(:),y(:,:) !x(i), y(i,column)
!        character (len=*), intent(in) :: Dir, xlabel, ylabel
!        type(str), intent(in)::title(:)
!        character (len=255):: plot_title,command
!        character (len=*),optional :: plot_title_in
!        logical:: display,output
!        plot_title=title(1)%d
!        if (present(plot_title_in)) plot_title=plot_title_in
!        call gnuCommandFile(Dir, plot_title,title, xlabel, ylabel, display,output)
!        call gnuDataFile(Dir, x,y)
!        write(command,*)'cd '//trim(Dir)//' && gnuplot -persist gnuCF_f90.plt'//' &'
!        call execute_command_line(trim ( command ))
!        call gnuPlot_clean(Dir)
!    endsubroutine gnuPlot2D_xy_list
!
!    subroutine gnuPlot_clean(Dir)
!        character (len=*):: Dir
!        integer::io
!        open(unit=2, iostat=io,file=trim(Dir)//'/gnuCF_f90.plt', status='old')
!        close(2, status='delete')
!        open(unit=2, iostat=io, file=trim(Dir)//'/gnuDF_f90.dat', status='old')
!        close(2, status='delete')
!        !call execute_command_line('cd '//trim(Dir)//' && del gnuCF_f90.plt gnuDF_f90.dat'//' &')
!    end subroutine gnuPlot_clean
!
!    subroutine gnuCommandFile(Dir, plot_title,title, xlabel, ylabel,display,output)
!        type(str),intent(in)::title(:)
!        character (len=*),intent(in) :: Dir, plot_title, xlabel, ylabel
!        logical,intent(in):: display,output
!
!        open(10,file= trim(Dir)//'/gnuCF_f90.plt',status='replace')
!        write(10,*) 'set xlabel '//'"'//trim(xlabel)//'"'
!        write(10,*) 'set ylabel '//'"'//trim(ylabel)//'"'
!        if (.not. display .and. output) then
!            write(10,*) 'set terminal pdfcairo'
!            write(10,*) 'set output "' //trim(plot_title) //'.pdf"'
!            call plot_all_col
!            close(10)
!            return
!        endif
!        !write(10,*) 'set terminal wxt' !default terminal for Windows
!        call plot_all_col
!        if (output) then
!            write(10,*) 'pause -1'
!            write(10,*) 'set terminal pdfcairo'
!            write(10,*) 'set output "' //trim(title(1)%d) //'.pdf"'
!            write(10,*) 'replot'
!        endif
!        close(10)
!        contains
!        subroutine plot_all_col
!            integer::i
!            character (len=2)::ns
!            write(10,*) 'plot \'
!            do i=1,size(title)-1
!                write(ns,'(i0)')i+1
!                write(10,*) '"gnuDF_f90.dat" using 1:'//trim(ns)//' with lines title " '//trim(title(i)%d)//' "  ,\'
!            enddo
!            write(ns,'(i0)')size(title)+1
!            write(10,*) '"gnuDF_f90.dat" using 1:'//trim(ns)//' with lines title " '//trim(title(size(title))%d)//' " '
!        endsubroutine plot_all_col
!    endsubroutine gnuCommandFile
!
!    subroutine gnuDataFile(Dir, x,y)
!        real(8), intent(in):: x(:),y(:,:)
!        integer :: i
!        character (len=*) :: Dir
!        if (size(x) /= size(y,1)) pause "err: 2dPlot; array size mismatch"
!        open(05, file = trim(Dir)// '/gnuDF_f90.dat', status = 'replace')
!        do i = 1, size(x)
!            write(05,*) x(i), y(i,:)
!        enddo
!        close(05)
!    endsubroutine gnuDataFile
!endmodule gnuPlot_v1_0
