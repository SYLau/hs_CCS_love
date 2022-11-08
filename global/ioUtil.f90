!io Utilities
!io sample: write(*,fmt='("<i=",i0,", reals=",2(g0,1x),">")')100,100.0,100.d0

module io_mod
    implicit none
    private
    public::assert
    public::fileexist
    public::getUnit,writeTable,writeArray,reallocate_a,file_constTable,file_nCol,file_nRow,file_readTable,to_char
    public::append_a
    public::reverse_a
    public::createFile
    public::arth
    public::gen_list_ug_bdy,gen_list_ug_dt
    public::dp_seg_to_char,dp_exp_to_char


    interface writeTable
        module procedure writeTable_R, writeTable_RR
    endinterface writeTable

    interface writeArray
        module procedure writeArrayHeader,writeTitle,writeArray_1,writeArray_2,writeArray_3,writeArray_4 &
        ,writeArray_5,writeArray_6,writeArray_7,writeArray_8,writeArray_9, writeArray_10
    endinterface writeArray

    interface reallocate_a
        module procedure reallocate_iv_a,reallocate_rv_a, reallocate_rm_a,reallocate_cv_a, reallocate_cm_a
    endinterface reallocate_a

    interface append_a
        module procedure append_dv_a, append_cv_a
    end interface append_a

    interface reverse_a
        module procedure reverse_dv_a
    end interface reverse_a

    interface to_char
        module procedure dp_to_char, i_to_char
    end interface to_char

    interface arth
        module procedure arth_d
    end interface arth

    interface gen_list_ug_bdy
        module procedure gen_list_ug_bdy_dp
    end interface gen_list_ug_bdy

    interface gen_list_ug_dt
        module procedure gen_list_ug_dt_dp
    end interface gen_list_ug_dt

contains
    function getUnit()
        implicit none
        integer:: getUnit
        integer:: i, ios
        logical:: lopen

        getUnit = 0
        do i = 1, 999
            inquire(unit = i, opened = lopen, iostat = ios)
            if ( ios == 0 ) then
                if ( .not. lopen ) then
                    getUnit = i
                    return
                end if
            end if
        end do
        write(*,*) 'getUnit: no free units'
        stop
    endfunction getUnit

    subroutine assert(cond,msg)
        logical,intent(in)::cond
        character(len=*),intent(in)::msg
        if (.not.cond) call err_msg(msg)
    end subroutine assert

    subroutine err_msg(msg)
        character(len=*)::msg
        write(*,*) msg
        stop
    end subroutine err_msg

    !Modified from NRUtil
    function reallocate_iv_a(p,n)
        integer, dimension(:), allocatable :: p, reallocate_iv_a
        integer, intent(in) :: n
        integer :: nold,ierr
        allocate(reallocate_iv_a(n),stat=ierr)
        if (ierr /= 0) call err_msg('reallocate_iv_a: problem in attempt to allocate memory')
        if (.not. allocated(p)) return
        nold=size(p,1)
        reallocate_iv_a(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
    endfunction reallocate_iv_a
    function reallocate_rv_a(p,n)
        real(8), dimension(:), allocatable :: p, reallocate_rv_a
        integer, intent(in) :: n
        integer :: nold,ierr
        allocate(reallocate_rv_a(n),stat=ierr)
        if (ierr /= 0) call err_msg('reallocate_rv_a: problem in attempt to allocate memory')
        if (.not. allocated(p)) return
        nold=size(p,1)
        reallocate_rv_a(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
    endfunction reallocate_rv_a
    function reallocate_rm_a(p,n,m)
        real(8), dimension(:,:), allocatable :: p, reallocate_rm_a
        integer, intent(in) :: n,m
        integer :: nold,mold,ierr
        allocate(reallocate_rm_a(n,m),stat=ierr)
        if (ierr /= 0) call err_msg('reallocate_rm_a: problem in attempt to allocate memory')
        if (.not. allocated(p)) return
        nold=size(p,1)
        mold=size(p,2)
        reallocate_rm_a(1:min(nold,n),1:min(mold,m))=p(1:min(nold,n),1:min(mold,m))
        deallocate(p)
    endfunction reallocate_rm_a
    function reallocate_cv_a(p,n)
        complex(8), dimension(:), allocatable :: p, reallocate_cv_a
        integer, intent(in) :: n
        integer :: nold,ierr
        allocate(reallocate_cv_a(n),stat=ierr)
        if (ierr /= 0) call err_msg('reallocate_cv_a: problem in attempt to allocate memory')
        if (.not. allocated(p)) return
        nold=size(p,1)
        reallocate_cv_a(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
    endfunction reallocate_cv_a
    function reallocate_cm_a(p,n,m)
        complex(8), dimension(:,:), allocatable :: p, reallocate_cm_a
        integer, intent(in) :: n,m
        integer :: nold,mold,ierr
        allocate(reallocate_cm_a(n,m),stat=ierr)
        if (ierr /= 0) call err_msg('reallocate_cm_a: problem in attempt to allocate memory')
        if (.not. allocated(p)) return
        nold=size(p,1)
        mold=size(p,2)
        reallocate_cm_a(1:min(nold,n),1:min(mold,m))=p(1:min(nold,n),1:min(mold,m))
        deallocate(p)
    endfunction reallocate_cm_a

    !Append array
    function append_dv_a(a_old,a_add)
        real(8),dimension(:),intent(in)::a_old,a_add
        real(8),dimension(:),allocatable::append_dv_a
        integer :: nold,nadd,ierr
        nold=size(a_old,1)
        nadd=size(a_add,1)
        allocate(append_dv_a(nold+nadd),stat=ierr)
        if (ierr /= 0) then
            write(*,*)'append_dv_a: problem in attempt to allocate memory'
            stop
        endif
        append_dv_a(1:nold)=a_old
        append_dv_a(nold+1:nold+nadd)=a_add
    endfunction append_dv_a
    function append_cv_a(a_old,a_add)
        complex(8),dimension(:),intent(in)::a_old,a_add
        complex(8),dimension(:),allocatable::append_cv_a
        integer :: nold,nadd,ierr
        nold=size(a_old,1)
        nadd=size(a_add,1)
        allocate(append_cv_a(nold+nadd),stat=ierr)
        if (ierr /= 0) then
            write(*,*)'append_cv_a: problem in attempt to allocate memory'
            stop
        endif
        append_cv_a(1:nold)=a_old
        append_cv_a(nold+1:nold+nadd)=a_add
    endfunction append_cv_a

    !reverse array
    subroutine reverse_dv_a(arr)
        real(8),dimension(:),intent(inout)::arr
        real(8),dimension(size(arr))::temp
        integer::i,n
        n=size(arr)
        do i=1,n
            temp(i)=arr(n+1-i)
        end do
        arr=temp
    end subroutine reverse_dv_a

    !Read from data file
    subroutine file_nCol(path,nCol)
        implicit none
        character(len=*),intent(in)::path
        integer,intent(out)::nCol
        integer,parameter:: maxLineLen=1000,maxCol=30
        character(len=maxLineLen)::line, colTest(1:maxCol)
        integer::i, io, inFile
        inFile=getUnit()
        open(inFile,file=path,action='read',status='old')
        read(inFile,'(a)') line
        close(inFile)
        do i=1,maxCol
            read(line,*,iostat=io) colTest(1:i)
            if (io==-1) exit
        enddo
        nCol=i-1
    endsubroutine file_nCol

    subroutine file_nRow(path,nRow)
        implicit none
        character(len=*),intent(in)::path
        integer,intent(out)::nRow
        integer::io, inFile
        nRow=0
        inFile=getUnit()
        open(inFile,file=path,action='read',status='old')
        do
            read(inFile,*,iostat=io)
            if (io/=0) exit
            nRow=nRow+1
        enddo
        close(inFile)
    endsubroutine file_nRow

    subroutine file_readTable(path,table,skip) !read table with known nCol and nRow
        implicit none
        integer:: inFile
        character(len=*),intent(in):: path
        integer,optional,intent(in)::skip
        real(8), intent(out):: table(1:,1:)
        integer:: i,nCol,nRow,io
        nCol=size(table,1)
        nRow=size(table,2)
        inFile=getUnit()
        open(inFile,file=path,status='old',action='read')
        if (present(skip)) then
            do i=1,skip
                read(inFile,*,iostat=io)
                if (io/=0) call err_msg('err: iostat - file_readTable')
            enddo
        endif
        do i=1,nRow
            read(inFile,*,iostat=io) table(1:nCol, i)
            if (io/=0) call err_msg('err: iostat - file_readTable')
        enddo
        close(inFile)
    endsubroutine file_readTable

    subroutine file_constTable(path,table,skip) !read and allocate table according to nCol and nRow
        implicit none
        character(len=*),intent(in)::path
        real(8),allocatable,intent(inout)::table(:,:)
        integer,optional,intent(in)::skip
        integer::nCol,nRow
        if (allocated(table)) deallocate(table)
        call file_nCol(path,nCol)
        call file_nRow(path,nRow)
        if (present(skip))nRow=nRow-skip !title row & unit row
        allocate(table(1:nCol,1:nRow))
        call file_readTable(path,table,skip=skip)
    endsubroutine file_constTable

    !Writing data file
    subroutine writeTable_R(path, raw)
        implicit none
        character(len=*):: path
        real(8), intent(in):: raw(1:, 1:)
        integer:: outFile
        integer:: nCol, nRow, i
        outFile=getUnit()
        open(outFile, file= path,status = 'replace', action='write')
        nCol=size(raw,1)
        nRow=size(raw,2)
        do i= 1, nRow
            write(outFile, "(10(ES25.15))") raw(1:nCol,i)
        enddo
        close(outFile)
    endsubroutine writeTable_R

    subroutine writeTable_RR(path, raw1, raw2)
        implicit none
        character(len=*):: path
        real(8), intent(in):: raw1(1:), raw2(1:, 1:)
        integer:: outFile
        integer:: nCol, nRow, i

        outFile=getUnit()
        open(outFile, file= path,status = 'replace', action='write')
        nCol=size(raw2,1)
        nRow=size(raw1,1)
        do i= 1, nRow
            write(outFile, "(10(ES25.15))")raw1(i),  raw2(1:nCol,i)
        enddo
        close(outFile)
    endsubroutine writeTable_RR

    !Create file if not exist
    subroutine createFile(path,cc)
        character(len=*),intent(in)::path
        character(len=*),optional,intent(in):: cc
        integer::outFile
        logical::lexist
        integer:: i,pos,ncount

        outFile=getUnit()
        inquire(file=path,exist=lexist)
        if (.not. lexist) then
            open(outFile,file=path,status='replace')

            if (present(cc)) then
                ncount=0
                do i=1,len(cc)
                    if (cc(i:i)=='#') ncount=ncount+1
                end do
                if (ncount==0) then
                    write(outFile,*) cc
                    close(outFile)
                    return
                end if
                pos=0
                do i=1,ncount-1
                    pos=scan(cc(pos+1:),'#')+pos
                    write(outFile,'(a25)',advance='no') trim(cc(pos:scan(cc(pos+1:),'#')-1+pos))
                end do
                pos=scan(cc(pos+1:),'#')+pos
                write(outFile,'(a25)') trim(cc(pos:))
            end if

            close(outFile)
        endif
    end subroutine createFile

    function fileexist(path)
        character(len=*),intent(in)::path
        logical::fileexist
        inquire(file=path, exist=fileexist)
    end function fileexist

    !Write array
    subroutine writeArray_status(outFile,path,st,ac)
        integer,intent(in)::outFile
        character(len=*),intent(in)::path
        character(len=*),intent(in),optional::st,ac
        character(len=99)::status,access
        if (present(st)) then
            status=st
        else
            status='replace'
        end if
        if (present(ac)) then
            access=ac
        else
            access='sequential'
        end if
        open(outFile, file= path,status=status, action='write',access=access)
    end subroutine writeArray_status

    subroutine writeArrayHeader(path,cc,st,ac)
        implicit none
        character(len=*),intent(in):: path,cc
        character(len=*),intent(in),optional::st,ac
        integer:: outFile,i,pos,ncount
        outFile=getUnit()
        call writeArray_status(outFile,path,st,ac)
        ncount=0
        do i=1,len(cc)
            if (cc(i:i)=='#') ncount=ncount+1
        end do
        if (ncount==0) then
            write(outFile,*) cc
            close(outFile)
            return
        end if
        pos=0
        do i=1,ncount-1
            pos=scan(cc(pos+1:),'#')+pos
            write(outFile,'(a25)',advance='no') trim(cc(pos:scan(cc(pos+1:),'#')-1+pos))
        end do
        pos=scan(cc(pos+1:),'#')+pos
        write(outFile,'(a25)') trim(cc(pos:))
        close(outFile)
    endsubroutine writeArrayHeader

    subroutine writeTitle(path,cc,st,ac)
        implicit none
        character(len=*),intent(in):: path
        character(len=*),dimension(:),intent(in):: cc
        character(len=*),intent(in),optional::st,ac
        character(len=2),parameter::max_col='10'
        integer:: nRow, i,outFile
        call writeArray_status(outFile,path,st,ac)
        write(outFile,'('//trim(max_col)//'a25)') cc
        close(outFile)
    end subroutine writeTitle

    subroutine writeArray_1(path,a1,st,ac)
        implicit none
        character(len=*),intent(in):: path
        real(8), intent(in):: a1(:)
        character(len=*),intent(in),optional::st,ac
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, "(10(ES25.15))")a1(i)
        enddo
        close(outFile)
    endsubroutine writeArray_1
    subroutine writeArray_2(path,a1,a2,st,ac)
        implicit none
        character(len=*),intent(in):: path
        real(8), intent(in):: a1(:), a2(:)
        character(len=*),intent(in),optional::st,ac
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, "(10(ES25.15))")a1(i),a2(i)
        enddo
        close(outFile)
    endsubroutine writeArray_2
    subroutine writeArray_3(path,a1,a2,a3,st,ac)
        implicit none
        character(len=*),intent(in):: path
        real(8), intent(in):: a1(:),a2(:),a3(:)
        character(len=*),intent(in),optional::st,ac
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, "(10(ES25.15))")a1(i),a2(i),a3(i)
        enddo
        close(outFile)
    endsubroutine writeArray_3
    subroutine writeArray_4(path,a1,a2,a3,a4,st,ac)
        implicit none
        character(len=*),intent(in):: path
        real(8), intent(in):: a1(:),a2(:),a3(:),a4(:)
        character(len=*),intent(in),optional::st,ac
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, "(10(ES25.15))")a1(i),a2(i),a3(i),a4(i)
        enddo
        close(outFile)
    endsubroutine writeArray_4
    subroutine writeArray_5(path,a1,a2,a3,a4,a5,st,ac)
        implicit none
        character(len=*),intent(in):: path
        real(8), intent(in):: a1(:),a2(:),a3(:),a4(:),a5(:)
        character(len=*),intent(in),optional::st,ac
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, "(10(ES25.15))")a1(i),a2(i),a3(i),a4(i),a5(i)
        enddo
        close(outFile)
    endsubroutine writeArray_5
    subroutine writeArray_6(path,a1,a2,a3,a4,a5,a6,st,ac)
        implicit none
        character(len=*),intent(in):: path
        real(8), intent(in):: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:)
        character(len=*),intent(in),optional::st,ac
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, "(10(ES25.15))")a1(i),a2(i),a3(i),a4(i),a5(i),a6(i)
        enddo
        close(outFile)
    endsubroutine writeArray_6
    subroutine writeArray_7(path,a1,a2,a3,a4,a5,a6,a7,st,ac)
        implicit none
        character(len=*),intent(in):: path
        real(8), intent(in):: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:),a7(:)
        character(len=*),intent(in),optional::st,ac
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, "(10(ES25.15))")a1(i),a2(i),a3(i),a4(i),a5(i),a6(i),a7(i)
        enddo
        close(outFile)
    endsubroutine writeArray_7
    subroutine writeArray_8(path,a1,a2,a3,a4,a5,a6,a7,a8,st,ac)
        implicit none
        character(len=*),intent(in):: path
        real(8), intent(in):: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:),a7(:),a8(:)
        character(len=*),intent(in),optional::st,ac
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, "(10(ES25.15))")a1(i),a2(i),a3(i),a4(i),a5(i),a6(i),a7(i),a8(i)
        enddo
        close(outFile)
    endsubroutine writeArray_8
    subroutine writeArray_9(path,a1,a2,a3,a4,a5,a6,a7,a8,a9,st,ac)
        implicit none
        character(len=*),intent(in):: path
        real(8), intent(in):: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:),a7(:),a8(:),a9(:)
        character(len=*),intent(in),optional::st,ac
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, "(10(ES25.15))")a1(i),a2(i),a3(i),a4(i),a5(i),a6(i),a7(i),a8(i),a9(i)
        enddo
        close(outFile)
    endsubroutine writeArray_9

    subroutine writeArray_10(path,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,st,ac)
        implicit none
        character(len=*),intent(in):: path
        real(8), intent(in):: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:),a7(:),a8(:),a9(:), a10(:)
        character(len=*),intent(in),optional::st,ac
        integer:: nRow, i,outFile
        outFile=getUnit()
        !open(outFile, file= path,status = 'replace', action='write')
        call writeArray_status(outFile,path,st,ac)
        nRow=size(a1)
        do i= 1, nRow
            write(outFile, "(10(ES25.15))")a1(i),a2(i),a3(i),a4(i),a5(i),a6(i),a7(i),a8(i),a9(i),a10(i)
        enddo
        close(outFile)
    endsubroutine writeArray_10


    !Convert to characters
    function i_to_char(i_in)
        integer,intent(in)::i_in
        character(len=25)::i_to_char
        write(i_to_char,'(i25)')i_in
    end function i_to_char

    function dp_to_char(dp_in)
        real(8),intent(in)::dp_in
        character(len=25)::dp_to_char
        write(dp_to_char,'(es25.15)')dp_in
    end function dp_to_char

    function dp_seg_to_char(num,ni,nf)
        real(8),intent(in)::num
        integer,intent(in)::ni,nf
        character(len=nf-ni+1)::dp_seg_to_char
        integer::i,nth_dig(nf),last_dig
        integer::expo
        real(8)::mantissa

        if (ni<1 .or. nf>16) then
            write(*,*)'err: dp_seg_to_char real(8) no exceed length'
            stop
        end if
        expo=floor(log10(num))
        mantissa=num
        if (expo>=0) then
            do i=1,expo
                mantissa=mantissa/10
            end do
        else
            do i=1,-expo
                mantissa=mantissa*10
            end do
        end if
        do i=1,nf
            nth_dig(i)=floor(mantissa)
            mantissa=(mantissa-nth_dig(i))*10
        end do
        last_dig=floor(mantissa)
        if (last_dig>=5) then
            nth_dig(nf)=nth_dig(nf)+1
            do i=nf,2,-1
                if (nth_dig(i)==10) then
                    nth_dig(i)=0
                    nth_dig(i-1)=nth_dig(i-1)+1
                endif
            end do
        end if
        write(dp_seg_to_char,'(16(i1))')nth_dig(ni:nf)
    end function dp_seg_to_char

    function dp_exp_to_char(num)
        real(8),intent(in)::num
        character(len=25)::dp_exp_to_char
        !character(len=:),allocatable::dp_exp_to_char
        !allocate(character(digits(res))::dp_exp_to_char)
        write(dp_exp_to_char,'(i0)')floor(log10(num))
    end function dp_exp_to_char

    !Gen array

    function arth_d(first,increment,n) !From numerical recipe nrutil
	real(8), intent(in) :: first,increment
	integer, intent(in) :: n
	integer, parameter :: NPAR_ARTH=16,NPAR2_ARTH=8
	real(8), dimension(n) :: arth_d
	integer :: k,k2
	real(8) :: temp
	if (n > 0) arth_d(1)=first
	if (n <= NPAR_ARTH) then
		do k=2,n
			arth_d(k)=arth_d(k-1)+increment
		end do
	else
		do k=2,NPAR2_ARTH
			arth_d(k)=arth_d(k-1)+increment
		end do
		temp=increment*NPAR2_ARTH
		k=NPAR2_ARTH
		do
			if (k >= n) exit
			k2=k+k
			arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
			temp=temp+temp
			k=k2
		end do
	end if
	end function arth_d

    function gen_list_ug_bdy_dp(ti,tf,nTot)
        real(8),intent(in)::ti,tf
        integer,intent(in)::nTot
        real(8),allocatable::gen_list_ug_bdy_dp(:)
        integer::i
        allocate(gen_list_ug_bdy_dp(nTot))
        do i=1,nTot
            gen_list_ug_bdy_dp(i)=ti+i*(tf-ti)/nTot
        end do
    end function gen_list_ug_bdy_dp

    function gen_list_ug_dt_dp(ti,dt,nTot)
        real(8),intent(in)::ti,dt
        integer,intent(in)::nTot
        real(8),allocatable::gen_list_ug_dt_dp(:)
        integer::i
        allocate(gen_list_ug_dt_dp(nTot))
        do i=1,nTot
            gen_list_ug_dt_dp(i)=ti+i*dt/nTot
        end do
    end function gen_list_ug_dt_dp

endmodule io_mod
