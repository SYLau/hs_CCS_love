module nr_mod

    interface
        function golden(ax,bx,cx,func,tol,xmin)
        use nrtype
        real(dp), intent(in) :: ax,bx,cx,tol
        real(dp), intent(out) :: xmin
        real(dp) :: golden
        interface
            function func(x)
            use nrtype
            real(dp), intent(in) :: x
            real(dp) :: func
            end function func
        end interface
        end function golden
    end interface

    interface
        subroutine mnbrak(ax,bx,cx,fa,fb,fc,func)
            use nrtype
            real(dp), intent(inout) :: ax,bx
            real(dp), intent(out) :: cx,fa,fb,fc
            interface
                function func(x)
                    use nrtype
                    real(dp), intent(in) :: x
                    real(dp) :: func
                end function func
            end interface
        end subroutine mnbrak
    end interface

    interface
        subroutine polint(xa,ya,x,y,dy)
            use nrtype
            real(dp), dimension(:), intent(in) :: xa,ya
            real(dp), intent(in) :: x
            real(dp), intent(out) :: y,dy
        end subroutine polint
    end interface

    interface
        function rtsec(func,x1,x2,xacc)
            use nrtype
            real(dp), intent(in) :: x1,x2,xacc
            real(dp) :: rtsec
            interface
                function func(x)
                    use nrtype
                    real(dp), intent(in) :: x
                    real(dp) :: func
                end function func
            end interface
        end function rtsec
    end interface

    interface
        function rtbis(func,x1,x2,xacc)
        use nrtype; use nrutil, only : nrerror
        implicit none
        real(dp), intent(in) :: x1,x2,xacc
        real(dp) :: rtbis
        interface
            function func(x)
                use nrtype
                implicit none
                real(dp), intent(in) :: x
                real(dp) :: func
            end function func
        end interface
    end function rtbis
    end interface

    interface
        subroutine sort(arr)
        use nrtype
        real(dp), dimension(:), intent(inout) :: arr
        end subroutine sort
    end interface
endmodule nr_mod
